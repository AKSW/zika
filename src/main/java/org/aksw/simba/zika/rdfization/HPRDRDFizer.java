package org.aksw.simba.zika.rdfization;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.openrdf.model.Literal;
import org.openrdf.model.Model;
import org.openrdf.model.Statement;
import org.openrdf.model.URI;
import org.openrdf.model.ValueFactory;
import org.openrdf.model.impl.LinkedHashModel;
import org.openrdf.model.vocabulary.XMLSchema;
import org.openrdf.repository.Repository;
import org.openrdf.repository.RepositoryConnection;
import org.openrdf.repository.RepositoryException;
import org.openrdf.repository.RepositoryResult;
import org.openrdf.repository.sail.SailRepository;
import org.openrdf.rio.RDFFormat;
import org.openrdf.rio.RDFHandlerException;
import org.openrdf.rio.Rio;
import org.openrdf.sail.memory.MemoryStore;
import info.aduna.iteration.Iterations;
/**
 * RDFize ZIKA Virus HPRD input files into RDF 
 * @author Saleem
 *
 */
public class HPRDRDFizer {
	public static RepositoryConnection conn = null;
	public static String res = "http://insight-centre/zika/resource/";
	public static String vocab = "http://insight-centre/zika/schema/";
	public static Repository rep ;
	public static ValueFactory factory;

	public static void main(String[] args) throws RepositoryException, RDFHandlerException, IOException {
		String tissueFile = "D:/eclipse-mars-workspace/Zika/data/Step_ 5_HPRD_data/HPRD_TISSUE.txt";
		String ppiFile = "D:/eclipse-mars-workspace/Zika/data/Step_ 5_HPRD_data/HPRD_PPI.txt";
		rdfizeHPRD(ppiFile,tissueFile);
	}
	/**
	 * RDFize the HPRD PPI files 
	 * @param ppiFile PPI file
	 * @param tissueFile Tissue File
	 * @throws IOException
	 * @throws RepositoryException
	 * @throws RDFHandlerException
	 */
	public static void rdfizeHPRD(String ppiFile, String tissueFile) throws IOException, RepositoryException, RDFHandlerException {
		initializeRepo();
		System.out.println("parsing started...");
		String line;
		BufferedReader br = new BufferedReader(new FileReader(ppiFile));
		String header[] = br.readLine().split("\t");
		ConcurrentHashMap<String, HashSet<String>> hm = getMap(ppiFile);
		URI sbj, obj, pred;
		int counter = 0;
		for(String key:hm.keySet())
		{
			HashSet<String> setVals = hm.get(key);
			int count = 1;
			for(String setVal:setVals)
			{
				sbj = factory.createURI(res, key);
				pred = factory.createURI(vocab, "hasInteration");
				obj = factory.createURI(res,key+"-int"+count);
				String linePrts [] = setVal.split("\t");
				conn.add(sbj,pred,obj);
				conn.add(obj,factory.createURI(vocab, header[1].trim()), factory.createURI(res, "hprdID:"+linePrts[1].trim())); //interactor_1_hprd_id
				conn.add(obj,factory.createURI(vocab, header[2].trim()), factory.createURI(res, linePrts[2].trim())); //interactor_1_refseq_id
				conn.add(obj,factory.createURI(vocab, header[3].trim()), factory.createURI(res, "gene:"+linePrts[3].trim())); //interactor_2_geneSymbol
				conn.add(obj,factory.createURI(vocab, header[4].trim()), factory.createURI(res, "hprdID:"+linePrts[4].trim())); //interactor_2_hprd_id
				conn.add(obj,factory.createURI(vocab, header[5].trim()), factory.createURI(res, "uniProtID:"+linePrts[5].trim())); //interactor_2_refseq_id
				conn.add(obj,factory.createURI(vocab, header[6].trim()), factory.createLiteral(linePrts[6].trim(),XMLSchema.STRING)); //experiment_typeB
				conn.add(obj,factory.createURI(vocab, header[7].trim()), factory.createLiteral(linePrts[7].trim(),XMLSchema.STRING)); //refid
				attachTissueInfo(obj,linePrts[1].trim(),tissueFile) ;
				count++;

			}
			System.out.println(counter++);
		}
		writeModel();
	}
	/**
	 * Create Hashmap of file
	 * @param ppiFile PPI file
	 * @return Hashmap of file
	 * @throws IOException
	 */
	public static ConcurrentHashMap<String, HashSet<String>> getMap(String ppiFile) throws IOException {
		String line;
		BufferedReader br = new BufferedReader(new FileReader(ppiFile));
		ConcurrentHashMap<String, HashSet<String>> hm = new ConcurrentHashMap<String, HashSet<String>>();
		while ((line = br.readLine()) != null)
		{
			String [] linePrts = line.split("\t");
			String key = linePrts[0].trim();
			if(!hm.containsKey(key))
			{
				HashSet<String> setValues = new HashSet<String>();
				setValues.add(line);
				hm.put(key, setValues)	;
			}
			else
			{
				HashSet<String> setVal = hm.get(key);
				synchronized(setVal){
					setVal.add(line);
				}

			}
		}
		return hm;
	}
	/**
	 * Attach Tissue info into PPI file of HPRD
	 * @param obj Object
	 * @param hprdID HPRD Id to search
	 * @param tissueFile The target file to be searched
	 * @throws IOException
	 * @throws RepositoryException
	 */
	public static void attachTissueInfo(URI obj,String hprdID, String tissueFile) throws IOException, RepositoryException {

		BufferedReader br = new BufferedReader(new FileReader(tissueFile));
		String line;
		while ((line = br.readLine()) != null)
		{	
			if(line.contains(hprdID))
			{
				String [] linePrts = line.split("\t");
				conn.add(obj,factory.createURI(vocab, "expressionTerm"), factory.createLiteral(linePrts[3].trim(),XMLSchema.STRING)); 
				conn.add(obj,factory.createURI(vocab, "status"), factory.createLiteral(linePrts[4].trim(),XMLSchema.STRING)); 
				conn.add(obj,factory.createURI(vocab, "refID"), factory.createLiteral(linePrts[5].trim(),XMLSchema.STRING)); 
				break;
			}
		}
	}

	/**
	 * Initialize Sesame Repository
	 * @throws RepositoryException
	 */
	private static void initializeRepo() throws RepositoryException {
		rep = new SailRepository(new MemoryStore());
		rep.initialize();
		factory  = rep.getValueFactory();
		conn = rep.getConnection();

	}
	/**
	 * Write the model into file
	 * @throws RepositoryException
	 * @throws FileNotFoundException
	 * @throws RDFHandlerException
	 */
	public static void writeModel() throws RepositoryException, FileNotFoundException, RDFHandlerException {
		RepositoryResult<Statement> statements =  conn.getStatements(null, null, null, true);
		Model model = Iterations.addAll(statements, new LinkedHashModel());
		OutputStream output = new FileOutputStream(new File("zika-HPRD.ttl"));
		Rio.write(model, output, RDFFormat.TURTLE);
		System.out.println("RDFization is completed. Output written to zika-HPRD.ttl");

	}


}
