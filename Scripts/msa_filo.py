import subprocess
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator


def run_mafft_alignment_windows(input_fasta, output_msa, mafft_path):
    """
    Executa o MAFFT via subprocesso.
    """
    subprocess.run(
        [mafft_path, "--auto", input_fasta],
        stdout=open(output_msa, "w"),
        stderr=subprocess.DEVNULL,
        check=True
    )
    print(f"Alinhamento salvo em: {output_msa}")


def run_mafft_alignment_macOS(input_fasta, msa_fasta):
    """
    Executa o MAFFT via subprocesso (assume que o 'mafft' está instalado no sistema).
    """
    subprocess.run(
        ["mafft", "--auto", input_fasta],
        stdout=open(msa_fasta, "w"),
        stderr=subprocess.DEVNULL,
        check=True
    )

    print(f"Alinhamento salvo em: {msa_fasta}")

def convert_to_phylip(msa_fasta, phylip_out):
    """
    Converte FASTA alinhado para PHYLIP Relaxed.
    """
    alignment = AlignIO.read(msa_fasta, "fasta")
    AlignIO.write(alignment, phylip_out, "phylip-relaxed")
    print(f"Convertido para Phylip: {phylip_out}")


def generate_upgma_tree(phylip_file):
    """
    Gera e imprime a árvore UPGMA.
    """
    aln = AlignIO.read(phylip_file, "phylip-relaxed")

    aln_subset = aln[:10]

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aln_subset)

    constructor = DistanceTreeConstructor()
    upgmatree = constructor.upgma(dm)

    print("\n--- Árvore UPGMA (Objeto) ---")
    print(upgmatree)

    return upgmatree


def draw_tree_ascii(tree):
    """Desenha a árvore em ASCII."""
    print("\n--- Desenho ASCII ---")
    Phylo.draw_ascii(tree)

