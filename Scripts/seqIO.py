from Bio import Entrez
from Bio import SeqIO


def setup_entrez(email):
    """Configura o email para o Entrez."""
    Entrez.email = email


def fetch_and_parse_genome(accession_id):
    """
    Faz o download do genoma completo e imprime informações gerais.
    Retorna o record do GenBank.
    """
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession_id)
    gb_record = SeqIO.read(handle, "genbank")
    handle.close()

    # Lógica original de impressão (adaptada para ler do objeto gb_record)
    print(gb_record.id, gb_record.description[:100], "...")
    print("Sequence length: ", len(gb_record))
    print(len(gb_record.features), " features")
    print("from: ", gb_record.annotations["source"])

    return gb_record


def analyze_gene_features(gb_record):
    """
    Percorre as features do tipo 'gene' e imprime suas localizações e qualificadores.
    """
    for feature in gb_record.features:
        if feature.type == "gene":
            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            print("\nFEATURE: gene")
            print("  start:", start)
            print("  end:", end)

            for key, value in feature.qualifiers.items():
                print(f"  {key} :", value)


def fetch_specific_gene_region(accession_id, start, stop, output_gb, output_fasta):
    """
    Procura uma região específica do genoma e salva em GB e FASTA.
    """
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession_id,
        rettype="gb",
        retmode="text",
        seq_start=start,
        seq_stop=stop
    )

    gene_record = SeqIO.read(handle, "genbank")
    handle.close()

    SeqIO.write(gene_record, output_gb, "genbank")
    SeqIO.write(gene_record, output_fasta, "fasta")

    print(gene_record.id, gene_record.description[:100], "...")
    print("Sequence length:", len(gene_record.seq))
    print(len(gene_record.features), "features")
    print("from:", gene_record.annotations["source"])

    return gene_record


def fetch_protein(protein_id, output_fasta):
    """
    Procura uma proteína e salva em FASTA.
    """
    handle = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=protein_id)

    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()

    for prot_record in records:
        print(prot_record.id, prot_record.description[:100], "...")
        print("Sequence length:", len(prot_record.seq))
        print(len(prot_record.features), "features")
        print("from:", prot_record.annotations["source"])

        SeqIO.write(prot_record, output_fasta, "fasta")


def analyze_detailed_features(gene_record):
    """
    Imprime taxonomias, tipos de features e qualificadores detalhados
    """
    featgene = []
    for i in range(len(gene_record.features)):
        if gene_record.features[i].type == "locus_tag":
            featgene.append(i)
    print("number of gene annotated:", len(featgene))

    if "taxonomy" in gene_record.annotations:
        print(gene_record.annotations["taxonomy"])

    for feature in gene_record.features:
        print(feature.type)

    for feature in gene_record.features:
        print("FEATURE:", feature.type)
        for key, value in feature.qualifiers.items():
            print(" ", key, ":", value)
        print()

    for feature in gene_record.features:
        if "db_xref" in feature.qualifiers:
            print("FEATURE:", feature.type)
            for xref in feature.qualifiers["db_xref"]:
                print("  db_xref:", xref)
