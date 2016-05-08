package org.aksw.simba.zika.rdfization;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
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
 * RDFize ZIKA Virus PPI Clustering input files into RDF 
 * @author Saleem
 *
 */
public class ClusterRDFizer {
	public static RepositoryConnection conn = null;
	public static String res = "http://insight-centre/zika/resource/";
	public static String vocab = "http://insight-centre/zika/schema/";
	public static Repository rep ;
	public static ValueFactory factory;

	public static void main(String[] args) throws RepositoryException, RDFHandlerException, IOException {
		String clusterFile = "D:/eclipse-mars-workspace/Zika/data/Step_ 4_Clustering/protein_clusters.txt";
		String ppiFile = "D:/eclipse-mars-workspace/Zika/data/Step_ 4_Clustering/PPI.txt";
		rdfizeCluster(ppiFile,clusterFile);
	}
	/**
	 * RDFize the cluster and PPI files. Alos link them
	 * @param ppiFile PPI file
	 * @param clusterFile Cluster File
	 * @throws IOException
	 * @throws RepositoryException
	 * @throws RDFHandlerException
	 */
	public static void rdfizeCluster(String ppiFile, String clusterFile) throws IOException, RepositoryException, RDFHandlerException {
		initializeRepo();
		System.out.println("parsing started...");
		BufferedReader br = new BufferedReader(new FileReader(ppiFile));
		String line; String header[] = br.readLine().split("\t");
		int count = 0;
		line = br.readLine();
		String [] linePrts = line.split("\t");
		String intA = linePrts[0];
		String newIntA = intA;
		URI sbj, obj, pred;
		do
		{	
			linePrts = line.split("\t");
			newIntA = linePrts[0];

			if (intA.equals(newIntA))
				count++;
			else{
				count =1;
				intA = newIntA;}

			sbj = factory.createURI(res, intA.trim());
			pred = factory.createURI(vocab, "hasInteration");
			obj = factory.createURI(res,intA.trim()+"-int"+count);
			conn.add(sbj,pred,obj);
			conn.add(obj,factory.createURI(vocab, header[1].trim()), factory.createURI(res, "uniProtID:"+linePrts[1].trim())); //uniprotID_A
			conn.add(obj,factory.createURI(vocab, header[2].trim()), factory.createLiteral(linePrts[2].trim(),XMLSchema.STRING));//uniProtName_A
			conn.add(obj,factory.createURI(vocab, header[3].trim()), factory.createURI(res, "gene:"+linePrts[3].trim())); //geneName_A
			conn.add(obj,factory.createURI(vocab, header[4].trim()), factory.createURI(res, linePrts[4])); //interactor_B
			conn.add(obj,factory.createURI(vocab, header[5].trim()), factory.createURI(res, "uniProtID:"+linePrts[5].trim())); //uniprotID_B
			conn.add(obj,factory.createURI(vocab, header[6].trim()), factory.createLiteral(linePrts[6].trim(),XMLSchema.STRING)); //uniProtName_B
			conn.add(obj,factory.createURI(vocab, header[7].trim()), factory.createURI(res, "gene:"+linePrts[7].trim())); //geneName_B
			conn.add(obj,factory.createURI(vocab, header[8].trim()), factory.createLiteral(linePrts[8].trim(),XMLSchema.INTEGER)); //curationEvents	
			conn.add(obj,factory.createURI(vocab, header[9].trim()), factory.createLiteral(linePrts[9].trim(),XMLSchema.STRING)); //curationMethod
			conn.add(obj,factory.createURI(vocab, header[10].trim()), factory.createURI(res, "PMID:"+linePrts[10].trim())); //PMID
			conn.add(obj,factory.createURI(vocab, header[11].trim()), factory.createLiteral(linePrts[11]+linePrts[12]+linePrts[13],XMLSchema.STRING)); //publication
			attachClusterInfo(obj,intA,linePrts[4].trim(),linePrts[7].trim(),clusterFile) ;

		}while ((line = br.readLine()) != null);
		// System.out.println(count);
		writeModel();
	}
	/**
	 * Attach the cluster scores to PPI values
	 * @param obj Ojbect
	 * @param intA Ineraction A
	 * @param intB Interaction B
	 * @param geneB Gene B can be same of IntB
	 * @param clusterFile Cluster file to be searched
	 * @throws IOException
	 * @throws RepositoryException
	 */
	public static void attachClusterInfo(URI obj,String intA, String intB, String geneB, String clusterFile) throws IOException, RepositoryException {

		BufferedReader br = new BufferedReader(new FileReader(clusterFile));
		String line;
		while ((line = br.readLine()) != null)
		{	
			if(line.contains(intA) &&(line.contains(intB) || line.contains(geneB)))
			{
				String [] linePrts = line.split("\t");
				conn.add(obj,factory.createURI(vocab, "EBIName"), factory.createLiteral(linePrts[1].trim(),XMLSchema.INTEGER)); //ebi name
				conn.add(obj,factory.createURI(vocab, "clusterNo"), factory.createLiteral(linePrts[3].trim(),XMLSchema.INTEGER)); //cluster no
				conn.add(obj,factory.createURI(vocab, "SCPCScore"), factory.createLiteral(linePrts[4].trim(),XMLSchema.FLOAT)); 
				conn.add(obj,factory.createURI(vocab, "NumberOfNodes"), factory.createLiteral(linePrts[5].trim(),XMLSchema.INTEGER));
				conn.add(obj,factory.createURI(vocab, "algorithm"), factory.createLiteral(linePrts[6].trim(),XMLSchema.STRING)); 
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
		OutputStream output = new FileOutputStream(new File("zika-cluster.ttl"));
		Rio.write(model, output, RDFFormat.TURTLE);
		System.out.println("RDFization is completed. Output written to zika-cluster.ttl");

	}


}
