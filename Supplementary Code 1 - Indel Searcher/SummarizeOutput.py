import os, sys
import pandas as pd
import subprocess as sp
from pdb import set_trace

sOutput_dir = sys.argv[1]
sProject    = sys.argv[2]

def Parsing_summary():

    if not os.path.isdir("{outdir}/Summary_result".format(outdir=sOutput_dir)):
        os.mkdir("{outdir}/Summary_result".format(outdir=sOutput_dir))

    ## Partially Cat files into Final File
    nBins       = 20
    list_sFiles = os.listdir('%s/Summary' % sOutput_dir)
    list_sFiles = [sFile for sFile in list_sFiles if not sFile.endswith('.temp')]
    list_sFiles = [sFile for sFile in list_sFiles if not sFile.endswith('Summary_all.txt')]

    nTotalJobs = len(list_sFiles)
    list_nIndexKey = [[i + 1, int(nTotalJobs * (i) / nBins), int(nTotalJobs * (i + 1) / nBins)] for i in
                      range(nBins)]

    for nIndex, nFileS, nFileE in list_nIndexKey:
        list_sInFiles = list_sFiles[nFileS:nFileE]
        list_InFiles  = ['%s/Summary/%s' % (sOutput_dir, sFile) for sFile in list_sInFiles]

        print('Summary_Part%s %s' % (nIndex, len(list_InFiles)))
        sOutFile      = '%s/Summary/Summary_Part%s.temp' % (sOutput_dir, nIndex)

        sp.call('cat %s >  %s' % (' '.join(list_InFiles), sOutFile), shell=True )
    #loop END: nIndex, nFileS, nFileE
    sp.call('cat {outdir}/Summary/Summary_Part*.temp > {outdir}/Summary/Summary_all.txt'.format(outdir=sOutput_dir), shell=True)

    dfSummary                = pd.read_table('{outdir}/Summary/Summary_all.txt'.format(outdir=sOutput_dir), header=None)
    dfSummary.columns        = ['Barcode', 'Total', 'Insertion', 'Deletion', 'Complex']
    dfSummary                = dfSummary.groupby(['Barcode']).sum()
    dfSummary['Total_indel'] = dfSummary['Insertion'] + dfSummary['Deletion'] + dfSummary['Complex']
    dfSummary['IND/TOT']     = dfSummary['Total_indel'] / dfSummary['Total']
    dfSummary['IND/TOT'].fillna(0, inplace=True)
    dfSummary.to_csv('{outdir}/Summary_result/Summary_result.tsv'.format(outdir=sOutput_dir), sep='\t')


def Annotate_final_result():
	
    dfCount_INDEL = pd.read_table('{outdir}/result/freq/freq_result/Indel_summary.txt'.format(outdir=sOutput_dir), header=None)
    dfSummary     = pd.read_table('{outdir}/Summary_result/Summary_result.tsv'.format(outdir=sOutput_dir), index_col='Barcode')

    dfCount_INDEL[0] = dfCount_INDEL[0].str.replace('.INDEL_freq.txt', '')
    dfCount_INDEL.set_index(0, inplace=True)
    dfConcat_result = pd.concat([dfCount_INDEL, dfSummary.loc[:,['Total_indel', 'Total', 'IND/TOT']]],axis=1)
    dfConcat_result.dropna(inplace=True)
    dfConcat_result = dfConcat_result.reset_index()
    dfConcat_result = dfConcat_result.loc[:,['index','Total_indel', 'Total', 'IND/TOT', 1,2]]
    dfConcat_result.columns = ['Barcode', 'Total_indel', 'Total', 'IND/TOT', 'Match','Info']
    dfConcat_result = dfConcat_result.round(2)
    dfConcat_result.to_csv('{outdir}/Summary_result/{project}_Final_indel_result.tsv'.format(outdir=sOutput_dir, project=sProject),
                           sep='\t', index=False)

if __name__ == '__main__':
    Parsing_summary()
    Annotate_final_result()


