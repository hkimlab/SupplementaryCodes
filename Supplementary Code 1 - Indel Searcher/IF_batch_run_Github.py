import os, sys
import multiprocessing as mp
import subprocess as sp

##################
## User options ##
####################################################################################################################


nCORES      = 15     ## multi core number
iChunk_size = 400000 ## split files, must be multiples of 4. file size < 1G recommendation:40000, size > 1G recommendation:400000
sDATA_DIR   = '/scratch'
sINPUT_DIR  = '%s/Input'  % sDATA_DIR
sOUTPUT_DIR = '%s/Output' % sDATA_DIR
GRABIX      = <PATH to grabix>

class Path_info(object):

    def __init__(self, sProject, sLibrary, sFastqTag, sPamPos):
        self.sProject       = sProject
        self.sLibrary       = sLibrary
        self.sFastq_name    = sFastqTag
        self.sQual_cutoff   = '20'
        self.sGap_open      = '20'                     # needle aligner option
        self.sGap_extend    = '1'                      #
        self.sEnd_open      = '20'                     #
        self.sEnd_extend    = '1'                      #
        self.iInsertion_win = '4'                      # Insertion window 0,1,2,3,4
        self.iDeletion_win  = '4'                      # Deletion window 0,1,2,3,4
        self.sPAM_type      = 'Cas9'                   # CRISPR type : Cpf1(2 cleavages), Cas9(1 cleavage)
        self.sPAM_pos       = sPamPos                  # Barcode target position : Forward (barcode + target), Reverse (target + barcode)


        # immutable variable
        self.sBarcode       = '%s/Input/Reference/%s/Barcode.txt'            % (sDATA_DIR, self.sLibrary)
        self.sReference_seq = '%s/Input/Reference/%s/Reference_sequence.txt' % (sDATA_DIR, self.sLibrary)
        self.sTarget_seq    = '%s/Input/Reference/%s/Target_region.txt'      % (sDATA_DIR, self.sLibrary)
        self.sRef_path      = '%s/Input/Reference/%s/Reference.fa'           % (sDATA_DIR, self.sProject)  # reference e.g.Consensus.fa

        if not os.path.isdir('%s/Input/Reference/%s' % (sDATA_DIR, self.sProject)):
            os.mkdir('%s/Input/Reference/%s'         % (sDATA_DIR, self.sProject))

        self.sInput_file = '%s/Input/FASTQ/%s/%s.fastq' % (sDATA_DIR, self.sProject, self.sFastq_name)
        self.sInput_list = '%s/Input/FASTQ/%s/%s.txt'   % (sDATA_DIR, self.sProject, self.sFastq_name)  # splited input file names for fastq processing

        self.sOutput_dir = '%s/Output/%s' % (sDATA_DIR, self.sProject)
        if not os.path.isdir(self.sOutput_dir): os.mkdir(self.sOutput_dir)
        if not os.path.isdir(self.sOutput_dir+'/Summary'): os.mkdir(self.sOutput_dir+'/Summary')                                        
        if not os.path.isdir(self.sOutput_dir+'/result'): os.mkdir(self.sOutput_dir+'/result')
        if not os.path.isdir(self.sOutput_dir+'/result/freq'): os.mkdir(self.sOutput_dir+'/result/freq')
        if not os.path.isdir(self.sOutput_dir+'/result/freq/freq_result'): os.mkdir(self.sOutput_dir+'/result/freq/freq_result')

        self.sInput_path = os.path.dirname(self.sInput_file)
        self.sPair       = 'False'  # FASTQ pair: True, False
        self.sDir        = '/'.join(self.sInput_list.split('/')[:-1])


# for single node
class Single_node_controller(Path_info):

    def __init__(self, iCore, sProject, sLibrary, sFastqTag, sPamPos):
        super(Single_node_controller, self).__init__(sProject, sLibrary, sFastqTag, sPamPos)
        self.iCore = iCore

    def Split_file (self):

        sTempFile  = '%s/linecnt.txt' % self.sInput_path
        sScript    = 'wc -l %s > %s; ' % (self.sInput_file, sTempFile)
        os.system(sScript)

        iTotal_lines   = int([sReadLine.replace('\n', '').split(' ')[0] for sReadLine in open(sTempFile)][0])
        sFastqFile_zip = compress_and_index(self.sInput_file)  # Compressed Fastq File for Random Access


        list_nBins = [[i + 1, int(iTotal_lines * (i) / nCORES), int(iTotal_lines * (i + 1) / nCORES)-1] for i in
                      range(nCORES)]

        list_sFiles = []
        for nIndex, nLineS, nLineE in list_nBins:

            sFileName = '%s.fastq_%s.fq' % (self.sFastq_name, nIndex)

            print(sFileName, nLineS, nLineE)

            list_sFiles.append(sFileName)

            sOutFile = '%s/%s' % (self.sInput_path, sFileName)

            sScript = '%s grab %s %s %s > %s' % (GRABIX, sFastqFile_zip, nLineS, nLineE, sOutFile)

            os.system(sScript)
        # loop END: nIndex, nCntStop

        OutFile = open('%s/%s.txt' % (self.sInput_path, self.sFastq_name), 'w')
        for sFileName in list_sFiles:
            sOut = '%s\n' % sFileName
            OutFile.write(sOut)
        OutFile.close()
    #def END: Split_file_v2
    def Make_reference(self):
        with open(self.sBarcode) as Barcode, \
                open(self.sTarget_seq) as Target, \
                open(self.sReference_seq) as Ref, \
                open(self.sRef_path, 'w') as Output:

            lName = [sBar.replace('\n', '') + ':' + sTar for sBar, sTar in zip(Barcode, Target)]

            for i, sRow in enumerate(Ref):
                Output.write('>' + lName[i] + sRow)


    def Make_indel_searcher_CMD(self):
        lCmd = []
        sReverse = 'None'

        with open(self.sInput_list) as Input:
            for sFile in Input:
                lFile = sFile.replace('\n', '').split(' ')
                #print lFile, self.sDir
                sForward = self.sDir + '/' + lFile[0]
                if self.sPair == 'True':
                    sReverse = self.sDir + '/' + lFile[1]
 
                #sReverse = self.sDir + '/' + lFile[1]
                #print sForward, sReverse

                lCmd.append('./IF_searcher.py {forw} {reve} {ref} {pair} {GapO} {GapE} {EndO} {EndE} {Insertion_win} {Deletion_win} {PAM_type} {PAM_pos} {Qual} {outdir}'.format(forw=sForward, reve=sReverse,
                            ref=self.sRef_path, pair=self.sPair, GapO=self.sGap_open, GapE=self.sGap_extend, EndO=self.sEnd_open, EndE=self.sEnd_extend, Insertion_win=self.iInsertion_win, Deletion_win=self.iDeletion_win,
                            PAM_type=self.sPAM_type, PAM_pos=self.sPAM_pos, Qual=self.sQual_cutoff, outdir=self.sOutput_dir))
        return lCmd

    def Run_indel_freq_calculator(self):
        sp.call('./IF_calculator.py {outdir}'.format(outdir=self.sOutput_dir), shell=True)
        sp.call('./SummarizeOutput.py {outdir} {project}'.format(outdir=self.sOutput_dir, project=self.sProject), shell=True)

def compress_and_index (sInFile):
    sBGzipFile = '%s.bgz' % sInFile
    sScript    = 'cat %s | bgzip -c > %s;'  % (sInFile, sBGzipFile)
    sScript   += '%s index %s;'             % (GRABIX, sBGzipFile)
    os.system(sScript)
    return sBGzipFile
#def END: compress_and_index


def Run_indel_searcher(sCmd):
    sp.call(sCmd, shell=True)

def Run_multicore(lCmd, iCore):
    p = mp.Pool(iCore)
    p.map_async(Run_indel_searcher, lCmd).get()
    p.close()

def batch_run(sProject, sLibrary, sFastqTag):

    sPamPos     = 'Forward'
    Ins_single  = Single_node_controller(nCORES, sProject, sLibrary, sFastqTag, sPamPos)
    Ins_single.Split_file()
    Ins_single.Make_reference()

    lCmd        = Ins_single.Make_indel_searcher_CMD()
    Run_multicore(lCmd, nCORES)
    Ins_single.Run_indel_freq_calculator()
#def END: batch_run

def main():

    sProject  = 'Test'
    sLibrary  = 'Lib1'
    sFastqTaq = 'Test.fastq'

    batch_run(sProject, sLibrary, sFastqTaq)

#def END: main



if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__

