import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def _set_ArchR(genome: str, 
        seed: int = 2022) -> None:
    '''
    Basic settings for running ArchR
    ---
    Input: 
        - genome: mm9, mm10, hg19, hg38
        - threads: number of threads used to run ArchR -> set to be 1 otherwise will have some problem
    '''
    ArchR = importr('ArchR') ## import ArchR package in python
    ArchR_thread_func = robjects.r['addArchRThreads']
    ArchR_thread_func(threads=1)  ## force to set to 1
    seed_func = robjects.r['set.seed']
    seed_func(seed)
    ArchR_genome_func = robjects.r['addArchRGenome']
    ArchR_genome_func(genome)


def _run_ArchR(input_files: List[str], sample_names: List[str], 
        output_dir: str = 'ArchR_output',
        min_tss: int = 4, min_frags: int = 1000,
        add_tilemat = False, add_genescore = True):
    '''
    Run ArchR on input files
    ---
    Input: 
        - input_files: list of input files
        - sample_names: corresponding input sample names
        - min_tss: TSS enrichment score threshold to filter low quality cells
        - min_frags: number of mapped fragment threshold to filter low quality cells
    '''
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
 
    ArchR_create_func = robjects.r['createArrowFiles']
    ArchR_arrowfiles = ArchR_create_func(inputFiles=robjects.StrVector(input_files), 
            sampleNames=robjects.StrVector(sample_names),
            minTSS=min_tss, minFrags=min_frags, 
            addTileMat=add_tilemat, 
            addGeneScoreMat=add_genescore)
    ArchR_proj_func = robjects.r['ArchRProject']
    ArchR_proj = ArchR_proj_func(ArrowFiles=robjects.StrVector(ArchR_arrowfiles),
            outputDirectory=output_dir, copyArrows=True) ## copy everything there and delete all files generated
    ArchR_save_func = robjects.r['saveArchRProject']
    ArchR_save_func(ArchRProj=ArchR_proj, outputDirectory=output_dir, load=False)

    ArchR_genescore_func = robjects.r['getMatrixFromProject']
    ArchR_genescore_mat = ArchR_genescore_func(ArchR_proj, useMatrix="GeneScoreMatrix")







if __name__ == '__main__':
    data_dir = "/home/wma36/gpu/data/Mousebrain_scATACseq/dscATACseq"

    import os, glob
    suffix = '.fragments.tsv.gz'
    input_files =  glob.glob(data_dir+os.sep+'*'+suffix)
    sample_names = [os.path.basename(x).replace(suffix, '') for x in input_files]
