from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO


def run_blast(input_fasta, xml_output):
    """
    Submete uma sequência de proteína para análise BLASTP remota na base de dados 'nr'.

    Lê o ficheiro FASTA de entrada, envia a sequência para o serviço NCBIWWW e grava o resultado XML devolvido no caminho especificado.

    Args:
        input_fasta (str): O caminho para o ficheiro FASTA com a sequência de consulta.
        xml_output (str): O caminho onde o ficheiro XML de resultados será guardado.
    """
    protein = SeqIO.read(input_fasta, "fasta")

    result_handle = NCBIWWW.qblast("blastp", "nr", protein.format("fasta"))

    with open(xml_output, "w") as save_file:
        save_file.write(result_handle.read())

    result_handle.close()
    print(f"XML salvo em: {xml_output}")


def parse_blast_xml(xml_input):
    """
    Processa um ficheiro de saída BLAST em formato XML e apresenta os dados.

    Carrega os registos do ficheiro, exibe o número total de 'hits' e detalha as informações principais do primeiro
    alinhamento encontrado (id, descrição, comprimento e número de HSPs)

    Args:
        xml_input (str): O caminho para o ficheiro XML do BLAST.

    Returns:
        BlastRecord: O objeto contendo todos os resultados da análise.
    """
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
    Extrai e ordena os melhores alinhamentos de um resultado BLAST.

    Exclui resultados associados ao fago T7 (para evitar auto-correspondência),
    seleciona apenas alinhamentos com um E-value inferior ou igual a 1e-20
    e devolve os 10 melhores, ordenados pela significância estatística.

    Args:
        xml_input (str): O caminho para o ficheiro XML do BLAST.

    Returns:
        list: Uma lista de tuplos (accession, definition, evalue, sbjct).
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
    Filtra rigorosamente os melhores alinhamentos do BLAST.

    Exclui o fago T7 e também registos que contenham "vector", "synthetic" ou "fusion" na definição. Mantém apenas
    os 10 melhores resultados com E-value inferior a 1e-20.

    Args:
        xml_input (str): O caminho para o ficheiro XML do BLAST.

    Returns:
        list: Uma lista filtrada dos melhores 'hits'.
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
    Exporta a lista de alinhamentos selecionados para um ficheiro FASTA.

    Percorre a lista de 'hits' fornecida e grava cada sequência num ficheiro,
    formatando o cabeçalho com o accession number.

    Args:
        hits (list): Lista de tuplos contendo as informações dos alinhamentos.
        output_filename (str): O caminho para o ficheiro de saída.
    """
    with open(output_filename, "w") as f:
        for acc, desc, evalue, seq in hits:
            header = f"{acc}"
            f.write(f">{header}\n{seq}\n\n")


def merge_query_and_hits(query_fasta, hits_fasta, final_output):
    """
    Combina a sequência de consulta original com os 'hits' encontrados.

    Lê a sequência original, atribui-lhe o identificador fixo "NP_041960_1",
    insere-a no topo da lista de sequências encontradas e grava o conjunto
    final num novo ficheiro FASTA.

    Args:
        query_fasta (str): O ficheiro FASTA da sequência original.
        hits_fasta (str): O ficheiro FASTA com os 'hits' do BLAST.
        final_output (str): O caminho para o ficheiro combinado final.
    """
    records = list(SeqIO.parse(hits_fasta, "fasta"))
    query = SeqIO.read(query_fasta, "fasta")

    query.id = "NP_041960_1"
    query.name = "NP_041960_1"
    query.description = ""

    records.insert(0, query)

    SeqIO.write(records, final_output, "fasta")
    print(f"Arquivo final gerado: {final_output}")

