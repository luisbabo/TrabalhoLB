import subprocess
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator


def run_mafft_alignment_windows(input_fasta, output_msa, mafft_path):
    """
    Executa o alinhamento de sequências múltiplas utilizando o software MAFFT.

    Esta função invoca o executável do MAFFT através de um subprocesso,
    para ambientes onde é necessário fornecer o caminho completo do executável (Windows). Utiliza a opção '--auto'
    para selecionar automaticamente a melhor estratégia de alinhamento e redireciona o resultado para o ficheiro de
    saída especificado.

    Args:
        input_fasta (str): O caminho para o ficheiro FASTA de entrada.
        output_msa (str): O caminho onde o alinhamento final será guardado.
        mafft_path (str): O caminho completo para o executável do MAFFT.
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
    Realiza o alinhamento de sequências via MAFFT assumindo instalação no sistema.

    Executa o comando 'mafft' diretamente da linha de comandos, partindo do princípio que o software se encontra nas
    variáveis de ambiente (PATH), o que é habitual em sistemas macOS e Linux. Aplica a configuração automática e grava
    o resultado no ficheiro de destino.

    Args:
        input_fasta (str): O caminho para o ficheiro FASTA com as sequências.
        msa_fasta (str): O caminho para o ficheiro de saída do alinhamento.
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
    Converte o formato de alinhamento de FASTA para PHYLIP Relaxed.

    Utiliza o módulo AlignIO para ler o alinhamento múltiplo gerado e transcreve os dados para o formato PHYLIP
    (versão 'relaxed'), necessário para a compatibilidade com diversas ferramentas de análise filogenética.

    Args:
        msa_fasta (str): O caminho do ficheiro de alinhamento original (FASTA).
        phylip_out (str): O caminho onde o ficheiro PHYLIP será criado.
    """
    alignment = AlignIO.read(msa_fasta, "fasta")
    AlignIO.write(alignment, phylip_out, "phylip-relaxed")
    print(f"Convertido para Phylip: {phylip_out}")


def generate_upgma_tree(phylip_file):
    """
    Calcula a matriz de distâncias e constrói a árvore filogenética UPGMA.

    Carrega o alinhamento em formato PHYLIP, seleciona as primeiras 10 sequências, calcula a distância e gera a árvore
    utilizando o algoritmo UPGMA.

    Args:
        phylip_file (str): O caminho para o ficheiro de alinhamento PHYLIP.

    Returns:
        Tree: O objeto da árvore filogenética gerada.
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
    """
    Gera uma representação visual da árvore filogenética em formato ASCII.

    Utiliza as ferramentas de visualização do módulo Phylo para desenhar a árvore diretamente na consola, permitindo
    a visualização da estrutura de ramificação.

    Args:
        tree (Tree): O objeto da árvore filogenética a visualizar.
    """
    Phylo.draw_ascii(tree)