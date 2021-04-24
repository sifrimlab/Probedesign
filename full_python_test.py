from wimpy import chunks
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from csv import DictWriter
from pyfaidx import Fasta


correctlength_seq = []
GC_lower_bound: 43
GC_upper_bound: 63
Temp_lower_bound: 66
Temp_upper_bound: 76
kmer_length: 30
kmer_overlap: 20
probe_length: 135


def split_complement(transcript_file):
### Read in transcript files(genes of interest) -> transcripts in ensembleID form eg. ENSG00000010404
        ### Read in gencode reference sequence and take the reverse complement of the transcripts in transcript_file(genes of interest)
        genes = Fasta('gencode.v33.pc_transcripts.fixed.fa')
        complement_sequence = []
        with open(transcript_file,"r") as transcripts:
                transcripts = transcripts.readlines()
                transcripts = [x.rstrip() for x in transcripts]
        # with open("complement_sequence.fasta","w") as output:
        for transcript in transcripts:
            # Testing different reading methods
        #        output.append(genes[transcript].complement)
        #        output.append(genes.get.seq(transcript,rc=True))
        #    complement_sequence.append(genes[transcript].complement)
            complement_sequence.append(genes.get.seq(transcript,rc=True))

    split_complement_sequence = chunks(complement_sequence,chunk_size=30,overlap=20)
    return(split_complement_sequence)

def find_barcode(transcript, barcode_file):
    readout_selected = []
    with open(barcode_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            #print(row)
            readout_selected = []
            if row["id"] == transcript:
                for (key, value) in row.items():
                    if value == "1":
                        readout_selected.append(key)
                return(readout_selected)
        return(False)

if __name__ == "__main__":
    # TODO Transcript names are added to the barcoding file to get corresponding bits and RScodes

    # READ CONFIG
    fw_primer = Seq(GGGCCACGTCCTCAATCGAC)
    rv_primer = Seq(CCCTCGCCAAGGTTCGCTAG)
    with open("Readout_probes_information.csv","r") as readout_file, open("barcodes_merfish.csv","r") as barcode_file,open("encoding_probes.fasta","w") as output_encode:

        results = {}
        with open(readout_file, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                #print(row)
                results[row["Readout probe name"]] = row["Sequence"]

        for transcript in transcripts:
            readout_ids = find_barcode(transcript, barcode_file)


            for record in SeqIO.parse(split_complement_sequence,"fasta"):
                    if record.id.split("_")[0] == transcript:
                        s = random.sample(readout_ids,3)
                        readout_seq = [results[i] for i in s]
                        record.seq = fw_primer +"A" + Seq(readout_seq[0]) + "A" + record.seq + "A" + Seq(readout_seq[1]) + "A" + Seq(readout_seq[2]) + "A" + rv_primer
                        record.id =  record.id +"_" + "_".join(s)
                        print(record.id)
                        record.description = ""
                        #print(record.seq)
                        SeqIO.write(record, output_encode, "fasta")

                        return


def filter_probes(encoding_probes,GC_lower_bound=43,GC_upper_bound=63,templowerbound=66,tempupperbound=67,probe_length=135):
    filters= {}
    with open("filtered_probes.fasta","w") as output_filtered:
        with open("probe_stats.csv","w", newline='') as filter_file:
            fieldnames =["probe_id","probelength","GC","meltingtemp", "sequence"]
            csv_writer=csv.DictWriter(filter_file,fieldnames=fieldnames)
            csv_writer.writeheader()
            # for record in SeqIO.parse("encoding_probes.fasta","fasta"):
            for record in SeqIO.parse(encoding_probes,"fasta"):
                filters["probelength"] = len(record.seq)
                filters["probe_id"] = (record.id)
                filters["GC"] = (GC(record.seq))
                filters["meltingtemp"] = (mt.Tm_NN(record.seq))
                filters["sequence"] = (record.seq)
                csv_writer.writerow(filters)
                #print(filters)

                if (    len(record.seq) == probe_length
                and GC(record.seq) > GClowerbound
                and GC(record.seq) < GCupperbound
                and float(('%0.2f' % mt.Tm_NN(record.seq))) > templowerbound
                and float(('%0.2f' % mt.Tm_NN(record.seq))) < tempupperbound
                ):
                correctlength_seq.append(record)
                SeqIO.write(record, output_filtered, "fasta")
                return

### blasting the filtered probes over the internet
from Bio.Blast import NCBIWWW
record = SeqIO.read("filtered_probes.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
#We need to be a bit careful since we can use result_handle.read() to read the BLAST output only once
#calling result_handle.read() again returns an empty string. so we save it then close it then read it in again

with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

result_handle.close()

result_handle = open("my_blast.xml")

### blasting the filtered_probes locally
from Bio.Blast.Applications import NcbiblastxCommandline
help(NcbiblastxCommandline)

blastx_cline = NcbiblastxCommandline(query="opuntia.fasta", db="nr", evalue=0.001,
                                      outfmt=5, out="opuntia.xml")
blastx_cline
NcbiblastxCommandline(cmd='blastx', out='opuntia.xml', outfmt=5, query='opuntia.fasta',
db='nr', evalue=0.001)
print(blastx_cline)
    blastx -out opuntia.xml -outfmt 5 -query opuntia.fasta -db nr -evalue 0.001
stdout, stderr = blastx_cline()


### Checking the results
from Bio.Blast import NCBIXML
E_VALUE_THRESH = 1e-20
for record in NCBIXML.parse(open("results.xml")):
    if record.alignments:
        print("\n")
        print("query: %s" % record.query[:100])
        for align in record.alignments:
           for hsp in align.hsps:
              if hsp.expect < E_VALUE_THRESH:
                 print("match: %s " % align.title[:100])
