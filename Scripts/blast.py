from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO


def run_blast(input_fasta, xml_output):
    """Lê o fasta e executa o blastp contra nr."""
    protein = SeqIO.read(input_fasta, "fasta")

    result_handle = NCBIWWW.qblast("blastp", "nr", protein.format("fasta"))

    with open(xml_output, "w") as save_file:
        save_file.write(result_handle.read())

    result_handle.close()
    print(f"XML salvo em: {xml_output}")


def parse_blast_xml(xml_input):
    """Imprime estatísticas do primeiro alinhamento."""
    with open(xml_input) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    print("number hits: ", len(blast_record.alignments))
    if blast_record.alignments:
        first_alignment = blast_record.alignments[0]
        print("FIRST ALIGNMENT:")
        print("Accession: " + first_alignment.accession)
        print("Hit id: " + first_alignment.hit_id)
        print("Definition: " + first_alignment.hit_def)
        print("Alignment length: ", first_alignment.length)
        print("Number of HSPs: ", len(first_alignment.hsps))

    return blast_record

def best_hits_no_filter(xml_input):
    """
    Retorna lista de melhores hits.
    """
    blast_record = NCBIXML.read(open(xml_input))
    best_hits = []

    for alignment in blast_record.alignments:
        if "Escherichia phage T7" in alignment.hit_def:
            continue

        hsp = alignment.hsps[0]

        if hsp.expect <= 1e-20:
            best_hits.append((
                alignment.accession,
                alignment.hit_def,
                hsp.expect,
                hsp.sbjct
            ))

    best_hits = sorted(best_hits, key=lambda x: x[2])[:10]

    for hit in best_hits:
        print("Accession:", hit[0])
        print("Definition:", hit[1])
        print("E-value:", hit[2])
        print()

    return best_hits


def best_hits_filter(xml_input):
    """
    Filtra os hits removendo T7, vetores, sintéticos e fusões.
    Retorna lista de melhores hits.
    """
    blast_record = NCBIXML.read(open(xml_input))
    best_hits = []

    for alignment in blast_record.alignments:
        if "Escherichia phage T7" in alignment.hit_def:
            continue
        if "vector" in alignment.hit_def or "synthetic" in alignment.hit_def or "fusion" in alignment.hit_def:
            continue

        hsp = alignment.hsps[0]

        if hsp.expect <= 1e-20:
            best_hits.append((
                alignment.accession,
                alignment.hit_def,
                hsp.expect,
                hsp.sbjct
            ))

    best_hits = sorted(best_hits, key=lambda x: x[2])[:10]

    for hit in best_hits:
        print("Accession:", hit[0])
        print("Definition:", hit[1])
        print("E-value:", hit[2])
        print()

    return best_hits


def save_hits_fasta(hits, output_filename):
    """
    Salva os hits.
    """
    with open(output_filename, "w") as f:
        for acc, desc, evalue, seq in hits:
            header = f"{acc}"
            f.write(f">{header}\n{seq}\n\n")


def merge_query_and_hits(query_fasta, hits_fasta, final_output):
    """Adiciona a query original ao topo da lista de hits."""
    records = list(SeqIO.parse(hits_fasta, "fasta"))
    query = SeqIO.read(query_fasta, "fasta")

    query.id = "NP_041960_1"
    query.name = "NP_041960_1"
    query.description = ""

    records.insert(0, query)

    SeqIO.write(records, final_output, "fasta")
    print(f"Arquivo final gerado: {final_output}")

