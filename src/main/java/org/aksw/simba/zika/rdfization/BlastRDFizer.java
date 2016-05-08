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
 * RDFize ZIKA Virus BLAST input files into RDF 
 * @author Saleem
 *
 */
public class BlastRDFizer {
	public static RepositoryConnection conn = null;
	public static String res = "http://insight-centre.org/zika/resource/";
	public static String vocab = "http://insight-centre.org/zika/schema/";
	public static Repository rep ;
	public static ValueFactory factory;

	public static void main(String[] args) throws RepositoryException, RDFHandlerException, IOException {
		String blastScoreFile = "D:/eclipse-mars-workspace/Zika/data/step1_blast/tordfize/JEKANCHX01R-Alignment-Scores.txt";
		String blastAlignFile = "D:/eclipse-mars-workspace/Zika/data/step1_blast/tordfize/Alignments.txt";
		String blastFullSeqFile = "D:/eclipse-mars-workspace/Zika/data/step1_blast/tordfize/Full seq of BLASt Matched Proteins.txt";
		rdfizeBlast(blastScoreFile,blastAlignFile,blastFullSeqFile);
	}
	/**
	 * RDFize BLAST output files pertaining to ZIKA reference number
	 * @param blastScoreFile Blast file containing the scores of difference reference number
	 * @param blastAlignFile File containing alignments. Note we have have separated the score file and alignment file of BLAST output
	 * @param blastFullSeqFile BLAST full sequece file
	 * @throws IOException
	 * @throws RepositoryException
	 * @throws RDFHandlerException
	 */
	public static void rdfizeBlast(String blastScoreFile, String blastAlignFile, String blastFullSeqFile) throws IOException, RepositoryException, RDFHandlerException {
		initializeRepo();
		System.out.println("parsing started...");
		BufferedReader br = new BufferedReader(new FileReader(blastScoreFile));
		String line; String header[] = br.readLine().split("\t");
		String dataset, refNo, proteinType, virus, score, eval;
		URI sbj, obj, pred; Literal literalObj;
		int count =1 ;
		while ((line = br.readLine()) != null)
		{	
			String linePrts[] = line.split("\t");
			String firstTuplePrts[] = linePrts[0].split("\\|"); 
			//System.out.print(firstTuplePrts.length);
			dataset = firstTuplePrts[0];
			//System.out.println(dataset);
			refNo = firstTuplePrts[1];
			//System.out.println(refNo);
			String secondTuplePrts[] = linePrts[1].replace("]", "").split("\\["); 
			proteinType = secondTuplePrts[0];
			if(linePrts[1].contains("["))
				virus = secondTuplePrts[1];
			else
				virus ="NULL";
			score = linePrts[2];
			eval = linePrts[3];			
			sbj = factory.createURI(res, refNo);
			conn.add(sbj,factory.createURI(vocab, header[1]), factory.createLiteral(refNo,XMLSchema.STRING));
			conn.add(sbj,factory.createURI(vocab, header[0]), factory.createLiteral(dataset,XMLSchema.STRING));
			conn.add(sbj,factory.createURI(vocab, header[2]), factory.createLiteral(proteinType,XMLSchema.STRING));
			conn.add(sbj,factory.createURI(vocab, header[3]), factory.createLiteral(virus,XMLSchema.STRING));
			conn.add(sbj,factory.createURI(vocab, header[4]), factory.createLiteral(score,XMLSchema.LONG));
			conn.add(sbj,factory.createURI(vocab, header[5]), factory.createLiteral(eval,XMLSchema.FLOAT));
			addBlastAlignments(refNo, blastAlignFile);
			String fullSeq = getSbjFullSeq(refNo,blastFullSeqFile);
			conn.add(sbj,factory.createURI(vocab, "sbjFullSeq"), factory.createLiteral(fullSeq,XMLSchema.STRING));
			count++;
			System.out.println(count);
		}
		// System.out.println(count);
		writeModel();
	}
	/**
	 * Get Subject Full Subject Sequence for the given reference number
	 * @param refNo Reference no
	 * @param blastFullSeqFile The file to search
	 * @return Full Subject Sequence
	 * @throws IOException
	 */
	public static String getSbjFullSeq(String refNo, String blastFullSeqFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(blastFullSeqFile));
		String line;
		String fullSeq ="";
		//System.out.println(count);
		while ((line = br.readLine()) != null)
		{	
			if(line.contains(refNo))
			{
				while (!(line = br.readLine()).contains(">"))
				{
					fullSeq = fullSeq+line;
				}
				return fullSeq;

			}
		}
		return fullSeq;
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
		OutputStream output = new FileOutputStream(new File("zika-blast.ttl"));
		Rio.write(model, output, RDFFormat.TURTLE);
		System.out.println("RDFization is completed. Output written to zika-blast.ttl");

	}
	/**
	 * Add alignments into each reference number obtained from BLAST output file 
	 * @param refNo The reference number for which we need to add alignments
	 * @param blastAlignFile File containing alignments. Note we have have separated the score file and alignment file of BLAST output
	 * @throws IOException
	 * @throws RepositoryException
	 */
	public static void addBlastAlignments(String refNo, String blastAlignFile) throws IOException, RepositoryException {
		BufferedReader br = new BufferedReader(new FileReader(blastAlignFile));
		String line;
		//System.out.println(count);
		outerloop:
			while ((line = br.readLine()) != null)
			{	
				if(line.contains(refNo))
				{
					//System.out.println(refNo+" matched: "+ line);
					while ((line = br.readLine()) != null) 
					{
						if (line.contains("Identities")){
							String prts[] = line.split(",");
							URI sbj = factory.createURI(res, refNo);
							String[] predOjb = prts[0].split("=");
							conn.add(sbj,factory.createURI(vocab, predOjb[0].trim()), factory.createLiteral(predOjb[1].trim(),XMLSchema.STRING));
							predOjb = prts[1].split("=");
							conn.add(sbj,factory.createURI(vocab, predOjb[0].trim()), factory.createLiteral(predOjb[1].trim(),XMLSchema.STRING));
							predOjb = prts[2].split("=");
							conn.add(sbj,factory.createURI(vocab, predOjb[0].trim()), factory.createLiteral(predOjb[1].trim(),XMLSchema.STRING));
							br.readLine();
							while (!(line = br.readLine()).contains(">")) 
							{
								if(line.contains("Query")){
									prts = line.split("\t");
									String queryRef = refNo+"-"+prts[0]+":"+prts[1]+"-"+prts[3];
									URI object = factory.createURI(res, queryRef);
									conn.add(sbj,factory.createURI(vocab, "hasQueryRef"),object);
									conn.add(object,factory.createURI(vocab, "hasTargetSeq"),factory.createLiteral(prts[2].trim(),XMLSchema.STRING));
									br.readLine();
									line = br.readLine();
									prts = line.split("\t");
									String sbjRef = queryRef+"-"+prts[0]+":"+prts[1]+"-"+prts[3];
									URI nobject = factory.createURI(res, sbjRef);
									conn.add(object,factory.createURI(vocab, "hasSbjRef"),nobject);
									conn.add(nobject,factory.createURI(vocab, "hasTargetSeq"),factory.createLiteral(prts[2].trim(),XMLSchema.STRING));

								}

							}
							break outerloop;
						}


					}

				}
			}
	}

}
