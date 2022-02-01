import os, sys, glob, shutil
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
        add_tilemat = False, add_genescore = True,
        save_proj = True):
    '''
    Run ArchR on input files
    ---
    Input: 
        - input_files: list of input files
        - sample_names: corresponding input sample names
        - min_tss: TSS enrichment score threshold to filter low quality cells
        - min_frags: number of mapped fragment threshold to filter low quality cells
        - add_tilemat: add 500bps tile matrix to the object
        - add_genescore: add gene score matrix to the object
        - save_proj: save ArchR project as RDS

    '''
    if len(input_files) != len(sample_names):
        ## exit the program
        sys.exit("Length of input files and sample names are not equal. Please make sure these two arguments are corresponding to each other.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
 
    ArchR_create_func = robjects.r['createArrowFiles']
    ArchR_arrowfiles = ArchR_create_func(inputFiles=robjects.StrVector(input_files), 
            sampleNames=robjects.StrVector(sample_names),
            minTSS=min_tss, minFrags=min_frags, 
            addTileMat=add_tilemat, 
            addGeneScoreMat=add_genescore)
    if len(ArchR_arrowfiles) != len(input_files):
        ## TODO: add warning
        pass

    ArchR_proj_func = robjects.r['ArchRProject']
    ArchR_proj = ArchR_proj_func(ArrowFiles=robjects.StrVector(ArchR_arrowfiles),
            outputDirectory=output_dir, copyArrows=True) ## copy everything there and delete all files generated
    if save_proj:
        ArchR_save_func = robjects.r['saveArchRProject']
        ArchR_save_func(ArchRProj=ArchR_proj, outputDirectory=output_dir, load=False)

    ArchR_genescore_func = robjects.r['getMatrixFromProject']
    ArchR_genescore_mat = ArchR_genescore_func(ArchR_proj, useMatrix="GeneScoreMatrix")

    mm, base = importr('Matrix'), importr("base") ## import Matrix package in python
    dollar = base.__dict__["$"]
    mm_write_func, write_func = robjects.r['writeMM'], robjects.r['write']
    assay_func, rowData_func, colnames_func = robjects.r['assays'], robjects.r['rowData'], robjects.r['colnames']
    mm_write_func(dollar(assay_func(ArchR_genescore_mat), 'GeneScoreMatrix'),
            output_dir+os.sep+'ArchR_genescore.mtx')
    write_func(dollar(rowData_func(ArchR_genescore_mat), 'name'),
            output_dir+os.sep+'ArchR_genescore_genes.tsv')
    write_func(colnames_func(ArchR_genescore_mat),
            output_dir+os.sep+'ArchR_genescore_barcodes.tsv')

    ## deal with files
    cur_dir = os.getcwd()
    qc_dir = cur_dir+os.sep+'QualityControl'
    if os.path.exists(qc_dir):
        shutil.move(qc_dir, output_dir+os.sep+'QualityControl') ## check filter
    log_dir = cur_dir+os.sep+'ArchRLogs'
    if os.path.exists(log_dir):
        shutil.move(log_dir, output_dir+os.sep+'ArchRLogs') ## check log for ArchR
    tmp_dir = cur_dir+os.sep+'tmp'
    if os.path.exists(tmp_dir) and len(os.listdir(tmp_dir)) == 0:
        os.rmdir(tmp_dir)  ## remove tmp directories
    else:
        ## TODO: add warining
        pass

    arrowfiles_dir = output_dir+os.sep+'ArrowFiles'
    if os.path.exists(arrowfiles_dir) and len(os.listdir(arrowfiles_dir)) != 0:
        cur_arrowfiles = glob.glob(cur_dir+os.sep+'*.arrow')
        for x in cur_arrowfiles: os.remove(x)  ## remove arrow files
    else:
        ## TODO: add warning
        pass

if __name__ == '__main__':
    data_dir = "/home/wma36/gpu/data/Mousebrain_scATACseq/dscATACseq"

    suffix = '.fragments.tsv.gz'
    input_files =  glob.glob(data_dir+os.sep+'*'+suffix)
    sample_names = [os.path.basename(x).replace(suffix, '') for x in input_files]

    _set_ArchR(genome='mm10', seed=2022)
    _run_ArchR(input_files, sample_names, output_dir=data_dir)
