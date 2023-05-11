#!/home/jinmanlab/bin/python3

import os, sys, random, pickle, time, subprocess, json, re
import subprocess as sp
import multiprocessing as mp
import numpy as np
from pdb import set_trace
from datetime import datetime

'''
< Method >
forward : 20bp + NGG
reverse : CCN + 20bp
N -> all base is available.

reference:

hg19 sequence have 'N'.
So if you use it, should filter it before calculation.
'''

sBASE_DIR   = '/data/scripts/gof_crisprcancer'
sINPUT_DIR  = '%s/Input'        % sBASE_DIR
sOUTPUT_DIR = '%s/Output'       % sBASE_DIR
sREF_DIR   = '/data/reference_genome'

sBE_DIR     = '%s/BaseEditing'  % sOUTPUT_DIR
sTEMP_DIR   = '%s/temp'         % sOUTPUT_DIR

sDB_DIR     = '/data/scripts/bin/annovar/humandb'
sAVINPUT    = '/data/scripts/bin/annovar/convert2annovar.pl'
#sANNOVAR    = '/data/scripts/bin/annovar/table_annovar.pl'
sANNOVAR    = '/data/scripts/bin/annovar/annotate_variation.pl'

sTIME_STAMP = '%s'              % (time.ctime().replace(' ', '-').replace(':', '_'))

################################### MAIN ANALYSIS CONDITIONS #############################################
# BE:  C -> T (G -> A)
# ABE: A -> G (T -> C)
sANALYSIS    = 'ABE'
nWIN_SIZE    = 7  # 7 = 4-7 | 8 = 4-8 | 9 = 4-9
##########################################################################################################

os.makedirs(sINPUT_DIR,  exist_ok=True)
os.makedirs(sOUTPUT_DIR, exist_ok=True)
os.makedirs(sBE_DIR,     exist_ok=True)
os.makedirs(sTEMP_DIR,   exist_ok=True)


dict_sKEY       = {'BE':'C>T', 'ABE':'A>G'}
list_sKEYCHECK  = ['C>T', 'G>A', 'A>G', 'T>C']
list_sCHRIDs    = ['chr' + str(x) for x in range(1, 23)] + ['chrX']
list_nCHRIDs    = [str(x) for x in range(1, 23)] + ['X', 'Y', 'M']

dict_sALTKEYS   = {'BE': {'+':['C','T'], '-':['G', 'A']},
                   'ABE':{'+':['A','G'], '-':['T', 'C']}}


##region Paths - OLD
class Path():

    sPrevNGG_Dir                    = '%s/fromHKK'                          % sINPUT_DIR
    sNGG_output                     = '%s/NGG_result.txt'                   % sPrevNGG_Dir
    sExtract_sg_only_path           = '%s/NGG_sgRNA_only_result.txt'        % sPrevNGG_Dir

    # Inputs
    sCCDS                           = '%s/UCSC_CCDS_exon_hg19_merged.bed'   % sINPUT_DIR
    sClinvar                        = '%s/clinvar_GRCh37_20170905.vcf'      % sINPUT_DIR
    sGWAS                           = '%s/GWAS_UCSC_hg19.txt'               % sINPUT_DIR
    sEXAC                           = '%s/hg19_exac03.txt'                  % sINPUT_DIR

    #sOutput_result                = './Results/BaseEditing/NGG_with_CCDS_{strand}_result.txt'

    sSNP_overlap_result             = '%s/SNP_{strand}_overlap_CCDS_result.txt'      % sBE_DIR
    sSNP_multi_overlap_result       = '%s/SNP_{strand}_multi_overlap_CCDS_result.txt'% sBE_DIR

    # Other Trials by Jaewoo
    sExtended_result                = '%s/NGG_with_extended.txt'                     % sBE_DIR
    sGWAS_overlap_result            = '%s/GWAS_{strand}_overlap_CCDS_result.txt     '% sBE_DIR


    '''NGG_output
    chrX    60005   60027   NNNNNNCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC     R
    chrX    60006   60028   NNNNNCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT     R
    chrX    60039   60061   CCCTAACCCTCTGAAAGTGGACCTATCAGCAGGATGTGGGTGG     F
    chrX    60045   60067   CCCTCTGAAAGTGGACCTATCAGCAGGATGTGGGTGGGAGCAG     F
    '''
##endregion

##region Colors

class cColors:
   PURPLE   = '\033[95m'
   CYAN     = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE     = '\033[94m'
   GREEN    = '\033[92m'
   YELLOW   = '\033[93m'
   RED      = '\033[91m'
   BOLD     = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'



def red(string, e=0):     return '\033[%s31m%s\033[m'%('' if e == 0 else '1;', string)
def green(string, e=0):   return '\033[%s32m%s\033[m'%('' if e == 0 else '1;', string)
def yellow(string, e=0):  return '\033[%s33m%s\033[m'%('' if e == 0 else '1;', string)
def blue(string, e=0):    return '\033[%s34m%s\033[m'%('' if e == 0 else '1;', string)
def magenta(string, e=0): return '\033[%s35m%s\033[m'%('' if e == 0 else '1;', string)
def cyan(string, e=0):    return '\033[%s36m%s\033[m'%('' if e == 0 else '1;', string)
def white(string, e=0):   return '\033[%s37m%s\033[m'%('' if e == 0 else '1;', string)

def get_color (cMir, sResidue):
    if sResidue == cMir.sAltNuc: sResidue = red(sResidue,1)
    elif sResidue == cMir.sRefNuc: sResidue = green(sResidue,1)
    else: sResidue = blue(sResidue,1)
    return sResidue
#def END: get_color
##endregion

## region class cVCFData

class cVCFData:
    def __init__(self):
        self.sPatID         = ''
        self.sChrID         = ''
        self.nPos           = 0
        self.sDBSNP_ID      = ''
        self.sRefNuc        = ''
        self.sAltNuc        = ''
        self.fQual          = 0.0
        self.sFilter        = ''
        self.sInfo          = ''
        self.sFormat        = ''
        self.list_sMisc     = []

        #Extra
        self.nClusterID     = 0

    #def END: __int__


def cVCF_parse_vcf_files (sVCFFile):
    if not os.path.isfile(sVCFFile):
        sys.exit('File Not Found %s' % sVCFFile)

    list_sOutput = []
    InFile       = open(sVCFFile, 'r')

    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | 1       | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat

        if sReadLine.startswith('#'): continue  # SKIP Information Headers
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cVCF                = cVCFData()
        cVCF.sChrID         = 'chr%s' % list_sColumn[0]

        try: cVCF.nPos      = int(list_sColumn[1])
        except ValueError: continue

        cVCF.sDBSNP_ID      = list_sColumn[2]
        cVCF.sRefNuc        = list_sColumn[3]
        cVCF.sAltNuc        = list_sColumn[4]
        cVCF.fQual          = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
        cVCF.sFilter        = list_sColumn[6]
        cVCF.sInfo          = list_sColumn[7]

        dict_sInfo          = dict([sInfo.split('=') for sInfo in cVCF.sInfo.split(';') if len(sInfo.split('=')) == 2])


        try: cVCF.sAlleleFreq    = float(dict_sInfo['AF_raw'])
        except ValueError: cVCF.sAlleleFreq = np.mean([float(f) for f in  dict_sInfo['AF_raw'].split(',')])

        list_sOutput.append(cVCF)
    #loop END: sReadLine
    InFile.close()

    return list_sOutput
#def END: cVCF_parse_vcf_files


def parse_vcf_stdout2 (sStdOut):

    list_sOutput   = []
    for sReadLine in sStdOut:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
        sReadLine           = str(sReadLine, 'UTF-8')
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cVCF                = cVCFData()
        if list_sColumn[0] == 'MT': continue
        if list_sColumn[0].startswith('<GL00'): continue

        dict_sChrKey = {'X':'23', 'Y':'24'}
        cVCF.sChrID  = 'chr%s' % list_sColumn[0]

        if list_sColumn[0] in ['X', 'Y']:
            cVCF.nChrID  = int(dict_sChrKey[list_sColumn[0]])
        else:
            cVCF.nChrID  = int(list_sColumn[0])
        cVCF.nPos           = int(list_sColumn[1])
        cVCF.sDBSNP_ID      = list_sColumn[2]
        cVCF.sRefNuc        = list_sColumn[3]
        cVCF.sAltNuc        = list_sColumn[4]
        cVCF.fQual          = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
        cVCF.sFilter        = list_sColumn[6]
        cVCF.sInfo          = list_sColumn[7]
        dict_sInfo          = dict([sInfo.split('=') for sInfo in cVCF.sInfo.split(';') if len(sInfo.split('=')) == 2])

        try: cVCF.fAlleleFreq    = float(dict_sInfo['AF_raw'])
        except ValueError: cVCF.fAlleleFreq = np.mean([float(f) for f in  dict_sInfo['AF_raw'].split(',')])

        list_sOutput.append(cVCF)
    #loop END: sReadLine

    return list_sOutput
#def END: parse_vcf_stdout

## endregion

## region class cCosmic
class cCOSMIC:
    def __init__(self):
        self.sGeneName   = ''
        self.sAccID      = ''
        self.nCDSLen     = 0
        self.sHGCNID     = ''   # SKIP for now
        self.sSample     = ''   # SKIP for now
        self.sSampleID   = ''   # SKIP for now
        self.sTumorID    = ''   # SKIP for now
        self.sPriSite    = ''   # primary site  ex) pancreas
        self.sSiteSub1   = ''   # SKIP for now
        self.sSiteSub2   = ''   # SKIP for now
        self.sSiteSub3   = ''   # SKIP for now
        self.sPriHist    = ''   # primary histology
        self.sHistSub1   = ''   # SKIP for now
        self.sHistSub2   = ''   # SKIP for now
        self.sHistSub3   = ''   # SKIP for now
        self.bGenomeWide = ''   # ex) y or n
        self.sMutaID     = ''   # SKIP for now
        self.sAltType    = ''   # ex) c.35G>T
        self.sRef        = ''
        self.sAlt        = ''
        self.sAAType     = ''   # ex) p.G12V
        self.sMutaDescri = ''   # ex) Substitution - Missense
        self.sMutaZygo   = ''   # SKIP for now
        self.bLOH        = ''   # loss of heterzygosity ex) y or n
        self.sGRCh       = ''   # Genome Version
        self.sGenicPos   = ''   # 17:7673781-7673781
        self.nChrID      = ''   # 17  X = 24 Y = 25
        self.sChrID      = ''   # chr17 or chrX and chrY
        self.sPos        = ''   # 7673781   1-based
        self.sStrand     = ''
        self.bSNP        = ''   # ex) y and n
        self.sDelete     = ''   # ex) PATHOGENIC
    #def END : __init__

def cCos_parse_cosmic_consensus (sChrID_Match):

    sCosmicFile  = '%s/COSMIC/cosmic_mutations_hg19_updated.tsv' % sREF_DIR
    list_sOutput = []
    InFile       = open(sCosmicFile, 'r', encoding='utf-8', errors='ignore')

    dict_sTest   = {}

    for i, sReadLine in enumerate(InFile):

        if sReadLine.startswith('Gene'):continue

        list_sColumn = sReadLine.strip('\n').split('\t')

        '''
        if i == 0:
            list_sHeader = list_sColumn
        elif i == 1:
            for i,(a,b) in enumerate(zip(list_sHeader, list_sColumn)):
                print('%s,%s,%s' % (i,a,b))
        else: break
        '''

        cCos             = cCOSMIC()
        cCos.sGeneName   = list_sColumn[0].upper()
        cCos.sAccID      = list_sColumn[1]
        cCos.nCDSLen     = list_sColumn[2]
        cCos.sHGNCID     = list_sColumn[3]
        cCos.sSample     = list_sColumn[4]
        cCos.sSampleID   = list_sColumn[5]
        cCos.sTumorID    = list_sColumn[6]
        cCos.sPriSite    = list_sColumn[7]
        cCos.sSiteSub1   = list_sColumn[8]
        cCos.sSiteSub2   = list_sColumn[9]
        cCos.sSiteSub3   = list_sColumn[10]
        cCos.sPriHist    = list_sColumn[11]
        cCos.sHistSub1   = list_sColumn[12]
        cCos.sHistSub2   = list_sColumn[13]
        cCos.sHistSub3   = list_sColumn[14]
        cCos.bGenomeWide = True if list_sColumn[15] == 'y' else False

        cCos.sCOSID      = list_sColumn[16]
        cCos.sOldMutaID  = list_sColumn[17]
        cCos.sMutaID     = list_sColumn[18]
        cCos.sAltType    = list_sColumn[19]
        cCos.sAAType     = list_sColumn[20]
        cCos.sMutaDescri = list_sColumn[21]
        cCos.sMutaZygo   = list_sColumn[22]
        cCos.bLOH        = True if list_sColumn[23] == 'y' else False
        cCos.sGRCh       = list_sColumn[24]
        cCos.sGenicPos   = list_sColumn[25]
        if not list_sColumn[25]: continue # Skip those w/o position information

        cCos.nChrID      = list_sColumn[25].split(':')[0]

        if cCos.nChrID not in dict_sTest:
            dict_sTest[cCos.nChrID] = 0
        dict_sTest[cCos.nChrID] += 1

        if cCos.nChrID not in ['24','25']:
            cCos.sChrID      = 'chr%s' % cCos.nChrID
        else:
            dict_sChrKey = {'24':'chrX', '25':'chrY'}
            cCos.sChrID  = dict_sChrKey[cCos.nChrID]
        #if END

        ### Separate by ChrID ###
        if cCos.sChrID != sChrID_Match: continue
        #########################

        list_sPosCheck   = list(set(list_sColumn[25].split(':')[1].split('-')))
        if len(list_sPosCheck) > 1:
            cCos.sPos    = list_sPosCheck[0]
        else:
            cCos.sPos    = ''.join(list_sPosCheck)
        #if END:

        cCos.sStrand     = list_sColumn[26]
        cCos.bSNP        = True if list_sColumn[27] == 'y' else False
        sKey             = cCos.sAltType[-3:]

        #### KEY CHECK ####
        if sKey != dict_sKEY[sANALYSIS]:continue
        ###################

        cCos.sMutaStrand     = list_sColumn[28]
        cCos.sFATHMMPred     = list_sColumn[29] if list_sColumn[26] else 'NA'
        cCos.sFATHMMScore    = list_sColumn[30] if list_sColumn[27] else 'NA'
        cCos.sMSStatus       = list_sColumn[31] if list_sColumn[27] else 'NA'
        cCos.sPMID           = list_sColumn[32]
        cCos.sIDStudy        = list_sColumn[33]
        cCos.sSampleType     = list_sColumn[34]
        cCos.sTumorOrigin    = list_sColumn[35]
        cCos.sAge            = list_sColumn[36]
        cCos.sHGVSP          = list_sColumn[37]
        cCos.HGVSC           = list_sColumn[38]
        cCos.HGVSG           = list_sColumn[39]

        list_sOutput.append(cCos)
    #loop END: i, sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cCos_parse_cosmic_consensus : list_sOutput : Size = %d' % (len(list_sOutput)))

    sOutFile = '%s/list_cCosmic_updated_%s_%s.data' % (sTEMP_DIR, sChrID_Match, sANALYSIS)
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(list_sOutput, OutFile)
    OutFile.close()
#def END: cCos_parse_cosmic_consensus
#class END: cCosmic
## endregion

## region class cFasta
re_nonchr = re.compile('[^a-zA-Z]')
class cFasta:
    def __init__(self, sRefFile):

        #V-S Check
        if not os.path.isfile(sRefFile):
           sys.exit('(): File does not exist')

        self.InFile     = open(sRefFile, 'r')
        self.sChrIDList = []
        self.nChromLen  = []
        self.nSeekPos   = []
        self.nLen1      = []
        self.nLen2      = []

        #V-S Check
        if not os.path.isfile('%s.fai'%sRefFile):
           sys.exit('.fai file does not exist')

        InFile = open('%s.fai' % sRefFile, 'r')
        for sLine in InFile:
           list_sColumn = sLine.strip('\n').split() # Goes backwards, -1 skips the new line character

           self.sChrIDList.append  (list_sColumn[0])
           self.nChromLen.append   (int(list_sColumn[1]))
           self.nSeekPos.append    (int(list_sColumn[2]))
           self.nLen1.append       (int(list_sColumn[3]))
           self.nLen2.append       (int(list_sColumn[4]))
        #loop END: sLINE
        InFile.close()
        self.sType = []
    #def END: __init_

    def fetch(self, sChrom, nFrom = None, nTo = None, sStrand = '+'):
        assert sChrom in self.sChrIDList, sChrom
        nChrom = self.sChrIDList.index(sChrom)

        if nFrom == None: nFrom = 0
        if nTo   == None: nTo = self.nChromLen[nChrom]
        #if nTo >= self.nChromLen[nChrom]: nTo = self.nChromLen[nChrom]-1

        assert(0 <= nFrom) and(nFrom < nTo) and(nTo <= self.nChromLen[nChrom])

        nBlank = self.nLen2[nChrom] - self.nLen1[nChrom]

        nFrom  = int(nFrom +(nFrom / self.nLen1[nChrom]) * nBlank) # Start Fetch Position

        nTo    = int(nTo   +(nTo   / self.nLen1[nChrom]) * nBlank) # End Fetch Position

        self.InFile.seek(self.nSeekPos[nChrom] + nFrom)            # Get Sequence

        sFetchedSeq = re.sub(re_nonchr, '', self.InFile.read(nTo - nFrom))

        if   sStrand == '+':
            return sFetchedSeq

        elif sStrand == '-':
            return reverse_complement(sFetchedSeq)

        else:
            sys.exit('Error: invalid strand')
        #if END: sStrand
    #def END: fetch
#class END: Fasta
## endregion



## region Util Functions

def copy_temp_core_script (sWorkDir):

    os.makedirs('%s/temp' % sWorkDir, exist_ok=True)
    os.system('cp %s/B_K1000E_additional_anals.py %s/temp/tmp_script_%s.py'
              % (sWorkDir, sWorkDir, sTIME_STAMP))

    return '%s/temp/tmp_script_%s.py' % (sWorkDir, sTIME_STAMP)
#def END: copy_temp_core_script

def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '.':'.', '*':'*',
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq       = list(sSeq) # Turns the sequence in to a gigantic list
    list_sSeq       = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]
#def END: reverse_complement

## endregion


def Extract_NGG_20mer_seq():

    assert os.path.isfile(Path.sNGG_output) == False, 'File already exist, clean the file'

    for sKey in ['X', 'Y', 'M'] + list(range(1,23)):
        sKey = 'chr'+str(sKey)
        print(sKey)

        sFor_seq = ''
        lRev_seq = []

        with open('/media/hkim/SubHardDisk/Pipeline/DNaseq_pipeline/Prerequisite/Genome/hg19/Splited/%s.fa' % sKey) as Chr:
            Chr.readline() # >chr1 : header skip.
            sFor_seq = Chr.read().upper().replace('\n', '')  # caution 'N'

        iFor_seq_len = len(sFor_seq)

        with open(Path.sNGG_output, 'a') as Output:

            for i in range(iFor_seq_len):
                sFor_target_seq = sFor_seq[i:i+23]
                sStart = i + 1  # convert one base
                sEnd   = i + 23

                if sFor_target_seq.find('N') == -1 and sFor_target_seq[-2:] == 'GG':    # This 'N' is in the hg19 reference.

                    '''
                    |--20bp--|
                    **********    N G G
                            -4   -3-2-1
                    |-GCcount-|
                    Until 4 decimal place. e.g. 0.0001
                    '''
                    sFor_target_seq_bothside_10bp = sFor_seq[i-10:i+33]

                    Output.write('\t'.join(map(str, [sKey, sStart, sEnd, sFor_target_seq_bothside_10bp, 'F']))+'\n')

                if sFor_target_seq.find('N') == -1 and sFor_target_seq[:2] == 'CC':

                    '''
                          |-20bp-|
                    C C N *********
                    0 1 2 3
                    Until 4 decimal place. e.g. 0.0001
                    '''
                    sRev_target_seq_bothside_10bp = sFor_seq[i-10:i+33]
                    Output.write('\t'.join(map(str, [sKey, sStart, sEnd, sRev_target_seq_bothside_10bp, 'R'])) + '\n')
            #end: for i
        #end: with open
    #end: for sKey
#def END: Extract_NGG_20mer_seq


def Extract_sgRNA_only():


    with open(Path.sNGG_output) as Input,\
        open(Path.sExtract_sg_only, 'w') as Output:

        for sRow in Input:
            lCol = sRow.replace('\n', '').split('\t')
            ssgRNA_additioanl_20bp = lCol[3]
            Output.write(ssgRNA_additioanl_20bp[10:-10]+'\n')
#def END: Extract_sgRNA_only


def set_intervalTree(sFile_path):

    '''
    For efficient searching algorithm, I use interval tree.
    dChr_Interval_CDS is consist of
    {'chr1': intervaltree object, 'chr2': interval...}
    intervaltree object : [geneSymbol, transcriptID, chrom, strand]
    intervaltree range  : [ExonStart:ExonEnd] This range is only on the CDS.
    '''
    dChr_Interval_CDS = {}

    for sKey in ['X', 'Y', 'M'] + range(1,23):
        sKey = 'chr' + str(sKey)
        dChr_Interval_CDS[sKey] = IntervalTree()

    ## Refflat_to_only_CDS_exon.py first
    with open(sFile_path) as Input:

        '''
        This is for intron PAM
         ------|   exon   |--------
                <-------NGG:F-->
          <--NGG:R------->       
        '''

        for sRow in Input:
            lCol   = sRow.replace('\n', '').split('\t')
            sChr   = lCol[0]
            iStart = int(lCol[1])
            iEnd   = int(lCol[2])
            if iEnd - iStart < 22:
                continue

            dChr_Interval_CDS[sChr][iStart :iEnd] = lCol  ## bed format.

    return dChr_Interval_CDS
#def END: set_intervalTree


def Set_interval_tree(sNGG_strand):  ## NGG with CDS, forward, reverse

    dChr_Interval_NGG_CDS = {}

    for sChrID in list_sCHRIDs:
        dChr_Interval_NGG_CDS[sChrID] = IntervalTree()

    sInFile  = '%s/NGG_with_CCDS_%s_result.txt' % (sBE_DIR, sNGG_strand)
    #if END:

    ## Refflat_to_only_CDS_exon.py first
    with open(sInFile) as Input:

        # Almost column is same. only last columns is different. The last columns is exon_start-exon_end on CDS.
        for sRow in Input:
            lCol     = sRow.replace('\n', '').split('\t')
            sChr     = lCol[0]
            iStart   = int(lCol[1]) - 1   ## NGG file 1base
            iEnd     = int(lCol[2])
            sStrand  = lCol[4]

            if sNGG_strand == 'F' and sStrand == 'F':
                dChr_Interval_NGG_CDS[sChr][iStart+3:iEnd-15] = lCol + [] ## bed format.
            elif sNGG_strand == 'R' and sStrand == 'R':
                dChr_Interval_NGG_CDS[sChr][iStart+15:iEnd-3] = lCol + []
            #if END:
        #loop END: sRow
    #with END:

    sOutFile = '%s/dChr_Interval_NGG_CDS_%s.data' % (sTEMP_DIR, sNGG_strand)
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dChr_Interval_NGG_CDS, OutFile)
    OutFile.close()
#def END: Set_interval_tree


def Set_interval_tree_v2(list_sParameters):  ## NGG with CDS, forward, reverse

    sChrID, sStrand = list_sParameters
    sOutObject      = IntervalTree()

    sInFile         = '%s/NGG_with_CCDS_%s_result.txt' % (sBE_DIR, sStrand)
    InFile          = open(sInFile, 'r')

    # Almost column is same. only last columns is different. The last columns is exon_start-exon_end on CDS.
    for sRow in InFile:

        lCol     = sRow.replace('\n', '').split('\t')
        sChr     = lCol[0]

        if sChrID != sChr: continue # Separate by match ChrIDs

        iStart   = int(lCol[1]) - 1   ## NGG file 1base
        iEnd     = int(lCol[2])
        sStrand  = lCol[4]

        if sStrand == 'F' and sStrand == 'F':
            sOutObject[iStart+3:iEnd-15] = lCol + [] ## bed format.
        elif sStrand == 'R' and sStrand == 'R':
            sOutObject[iStart+15:iEnd-3] = lCol + []
        #if END:
    #loop END: sRow
    InFile.close()

    sOutFile = '%s/dChr_Interval_NGG_CDS_%s_%s.data' % (sTEMP_DIR, sStrand, sChrID)
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(sOutObject, OutFile)
    OutFile.close()
#def END: Set_interval_tree_v2


def Write_CDS_NGG_data(objCheck, sRow, iStart, iEnd, CDS_only_output, CDS_without_output):

    if objCheck != []:
        iCDS_start = objCheck[0][0]
        iCDS_end = objCheck[0][1]

        iDist_st_st = iCDS_start - iStart
        iDist_st_en = iCDS_start - iEnd

        iDist_en_st = iCDS_end - iStart
        iDist_en_en = iCDS_end - iEnd
        CDS_only_output.write(sRow.replace('\n', '') + '\t' + '\t'.join(map(str, [iDist_st_st, iDist_st_en, iDist_en_st, iDist_en_en])) + '\n')

    elif objCheck == []:
        CDS_without_output.write(sRow)
#def END: Write_CDS_NGG_data


def filter_NGG_20mer(lParameter):

    dChr_Interval_CDS = lParameter[0]
    sNGG_Strand       = lParameter[1]

    print ('Filter NGG 20mer %s' % sNGG_Strand)

    sOut_wCDS  = '%s/NGG_with_CCDS_%s_result.txt'    % (sBE_DIR, sNGG_Strand)
    sOut_woCDS = '%s/NGG_without_CCDS_%s_result.txt' % (sBE_DIR, sNGG_Strand)

    with open(Path.sNGG_output) as NGG_input,\
        open(sOut_wCDS, 'w') as CDS_only_output,\
        open(sOut_woCDS, 'w') as CDS_without_output:

        for sRow in NGG_input:
            lCol   = sRow.replace('\n', '').split('\t')
            sChr   = lCol[0]
            iStart = int(lCol[1]) - 1   # one base position format, so -1
            iEnd   = int(lCol[2])
            sStrand = lCol[4]

            if sNGG_Strand == 'F' and sStrand == 'F':
                objCheck = sorted(dChr_Interval_CDS[sChr][iStart+3:iEnd-15])  # Check overlap exons on CDS.
                Write_CDS_NGG_data(objCheck, sRow, iStart, iEnd, CDS_only_output, CDS_without_output)

            elif sNGG_Strand == 'R' and sStrand == 'R':
                objCheck = sorted(dChr_Interval_CDS[sChr][iStart+15:iEnd-3])  # Check overlap exons on CDS.
                Write_CDS_NGG_data(objCheck, sRow, iStart, iEnd, CDS_only_output, CDS_without_output)
            #if END:
        #loop END:
    #with END:
#def END: filter_NGG_20mer


def Make_hg19_dict():

    if not os.path.isfile('/media/hkim/SubHardDisk/Pipeline/DNaseq_pipeline/Input/Pickle/hg19.pickle'):

        dHg19 = {}

        for sChr in range(1,23)+['X','Y','M']:
            sChr = 'chr'+str(sChr)
            with open('/media/hkim/SubHardDisk/Pipeline/Genome_data/Splited_genome/hg19/%s.fa' % sChr) as hg19_chr:
                hg19_chr.readline()
                dHg19[sChr] = hg19_chr.read().upper().replace('\n', '')
        pickle.dump(dHg19, open('/media/hkim/SubHardDisk/Pipeline/DNaseq_pipeline/Input/Pickle/hg19.pickle', 'wb'))

        return dHg19

    else:

        dHg19 = pickle.load(open('/media/hkim/SubHardDisk/Pipeline/DNaseq_pipeline/Input/Pickle/hg19.pickle'), 'rb')

        return dHg19
#def END: Make_hg19_dict


def Extend_200bp_with_NGG_20mer(dHg19, iRandom_selection_num):

    with open(Path.sExtended_result, 'w') as Result:
        for sRow in sp.check_output('shuf -n {ran_num} {NGG_CDS_result}'.format(ran_num=iRandom_selection_num, NGG_CDS_result=Path.sOutput_result), shell=True).split('\n'):
            # set_trace()
            lCol   = sRow.replace('\n', '').split('\t')
            sChr   = lCol[0]
            try:
                iStart = int(lCol[1]) -1 # 1base position format, not bed
                iEnd   = int(lCol[2])
            except IndexError:
                print (sRow)
                print (lCol)
                continue

            ## extend 200bp
            iStart = iStart - 200
            iEnd   = iEnd   + 200

            ## Temporary code, sgRNA seq is extended +-10bp. so remove it.
            lCol[3] = lCol[3][10:-10]

            sExtended_guide = dHg19[sChr][iStart :iEnd]
            Result.write('\t'.join(lCol+[sExtended_guide])+'\n')
#def END: Extend_200bp_with_NGG_20mer


def Make_EXAC_dict():

    dEXAC = {sChrID:{} for sChrID in list_sCHRIDs}

    with open(Path.sEXAC) as EXAC:

        # EXAC 1 base pos
        EXAC.readline()  # skip header

        for sRow in EXAC:
            lCol   = sRow.replace('\n', '').split('\t')
            sChr   = 'chr'+lCol[0]
            iStart = int(lCol[1])
            iEnd   = int(lCol[2])

            if iStart - iEnd == 0:

                dEXAC[sChr][iStart] = lCol[5]
        #END loop:
    #END with:

    sOutFile = '%s/dict_EXAC.data' % sTEMP_DIR
    OutFile  = open(sOutFile, 'wb')
    pickle.dump(dEXAC, OutFile)
    OutFile.close()
#def END: Make_EXAC_dict


def make_GNOMAD_by_chr (sChrID):

    sGnomadVCF  = '%s/gnomad.exomes.vcf.gz'                              % sINPUT_DIR
    sOutFile    = '%s/gnomAD_exome_vcf_bychr/gnomad.exomes.chr%s.vcf'    % (sINPUT_DIR, sChrID)
    sOutFile2   = '%s/gnomAD_exome_vcf_bychr/gnomad.exomes.chr%s.vcf.gz' % (sINPUT_DIR, sChrID)

    sScript     = 'tabix %s %s > %s;' % (sGnomadVCF, sChrID, sOutFile)
    sScript    += 'bgzip -c %s > %s;' % (sOutFile, sOutFile2)
    sScript    += 'tabix -p vcf %s;'  % sOutFile2

    #print(sScript)
    os.system(sScript)
#def END: make_GNOMAD_by_chr


def Make_SNP_overlap_result(dTarget_interval, sChr, iPos, lCol, SNP_overlap_output, sRef_base, sAlt_base,
                            dMulti_SNP_region, sStrand, dEXAC):

    objCheck = sorted(dTarget_interval[sChr][iPos])  # Check overlap exons on CDS.

    if objCheck != []:
        # one-base pos
        try:
            lInfo = lCol[7].split(';')
            lMain_info = ['NaN','NaN','NaN','NaN','NaN']
            for sInfo in lInfo:
                sMain_info = sInfo.split('=')[0]

                if sMain_info == 'CLNSRCID':
                    lMain_info[0] = sInfo

                elif sMain_info == 'CLNSIG':
                    lMain_info[1] = sInfo

                elif sMain_info == 'CLNDBN':
                    lMain_info[2] = sInfo

                elif sMain_info == 'CAF':
                    lMain_info[3] = sInfo

                elif sMain_info == 'COMMON':
                    lMain_info[4] = sInfo

            lNGG       = objCheck[0][2]
            sNGG_start = lNGG[1]
            sNGG_end   = lNGG[2]
            sNGG_seq   = lNGG[3][10:33]

            try:
                sEXAC_freq = dEXAC[sChr][int(lCol[1])]
            except KeyError:
                sEXAC_freq = 'Na'

            if sStrand == 'F': # one base calculation
                sForward_count_base = sNGG_seq[3:8]
                sResult = '\t'.join([sChr, sNGG_start, sNGG_end, sNGG_seq]) + \
                          '\t' + '\t'.join([sChr, lCol[1], str(int(lCol[1]) - int(sNGG_start) + 1), lCol[2], sRef_base, sAlt_base] + [sEXAC_freq]) + \
                          '\t' + '\t'.join(lMain_info+[str(sForward_count_base.count('A'))]) + '\n'

            elif sStrand == 'R':
                sReverse_count_base = sNGG_seq[-8:-3]
                sResult = '\t'.join([sChr, sNGG_start, sNGG_end, sNGG_seq]) + \
                          '\t' + '\t'.join([sChr, lCol[1], str(int(sNGG_end) - int(lCol[1]) + 1), lCol[2], sRef_base, sAlt_base] + [sEXAC_freq]) + \
                          '\t' + '\t'.join(lMain_info+[str(sReverse_count_base.count('T'))]) + '\n'

            SNP_overlap_output.write(sResult)
            lResult = sResult.replace('\n', '').split('\t')

            try:
                _ = dMulti_SNP_region[sChr + '_' + sNGG_start + '_' + sNGG_end]
                dMulti_SNP_region[sChr + '_' + sNGG_start + '_' + sNGG_end][0] += 1
                dMulti_SNP_region[sChr + '_' + sNGG_start + '_' + sNGG_end][1].append(lResult)
            except KeyError:
                dMulti_SNP_region[sChr + '_' + sNGG_start + '_' + sNGG_end] = [1, [lResult]]

            return lResult

        except IndexError:
            set_trace()
    #if END
#def END: Make_SNP_overlap_result


def Check_overlap_SNP(dTarget_interval, sStrand, dEXAC, sAnalysis):

    '''
    base editing target
    4,5,6,7base in NGG
    overlap snp.
    :param dTarget_interval:
    :return:
    '''

    with open(Path.sClinvar) as SNP,\
        open(Path.sSNP_overlap_result.format(strand=sStrand), 'w') as SNP_overlap_output,\
        open(Path.sSNP_multi_overlap_result.format(strand=sStrand), 'w') as SNP_multi_overlap_output:

        SNP_overlap_output.write('Chromosome\tNGG_start\tNGG_end\tNGG_seq\tClinvar_chr\tClinvar_pos\tBase_num\tDBSNP_ID\tRef\tAlt\tCLNSRCID\tCLNSIG\tCLNDBN\tCAF\tCOMMON\tMulti_base\n')
        SNP_multi_overlap_output.write('Chromosome\tNGG_start\tNGG_end\tNGG_seq\tClinvar_chr\tClinvar_pos\tBase_num\tDBSNP_ID\tRef\tAlt\tEXAC_ALL_freq\tCLNSRCID\tCLNSIG\tCLNDBN\tCAF\tCOMMON\tMulti_base\tMulti_SNP\n')
        dMulti_SNP_region = {}

        for sRow in SNP:                      # Clinvar 1base.
            if sRow[0] == '#': continue
            lCol = sRow.replace('\n','').split('\t')
            sChr = 'chr' + lCol[0]
            if sChr == 'chrMT': continue

            iPos = int(lCol[1]) - 1          # 1base -> 0 base

            sRef_base = lCol[3]
            sAlt_base = lCol[4]

            if sStrand == 'F' and sRef_base == 'G' and sAlt_base == 'A':
                lResult = Make_SNP_overlap_result(dTarget_interval, sChr, iPos, lCol, SNP_overlap_output, sRef_base, sAlt_base,
                                                  dMulti_SNP_region, sStrand, dEXAC)

            elif sStrand == 'R' and sRef_base == 'C' and sAlt_base == 'T':
                lResult = Make_SNP_overlap_result(dTarget_interval, sChr, iPos, lCol, SNP_overlap_output, sRef_base, sAlt_base,
                                                  dMulti_SNP_region, sStrand, dEXAC)

        for sKey, llValue in dMulti_SNP_region.items():
            iCount = llValue[0]
            for lValue in llValue[1]:
                SNP_multi_overlap_output.write('\t'.join(lValue) + '\t'+ str(iCount) + '\n')
#def END: Check_overlap_SNP


def Check_overlap_GWAS(sStrand):

    """
    base editing target

    4,5,6,7base in NGG
    overlap snp.
    :param dTarget_interval:
    :return:
    """
    with open(Path.sSNP_multi_overlap_result.format(strand=sStrand)) as SNP,\
        open(Path.sGWAS) as GWAS,\
        open(Path.sGWAS_overlap_result.format(strand=sStrand), 'w') as GWAS_anno_output:

        #SNP_multi_overlap_output.write('Chromosome\tNGG_start\tNGG_end\tNGG_seq\tClinvar_chr\tClinvar_pos\tBase_num\tDBSNP_ID\tRef\tAlt\tCLNSRCID\tCLNSIG\tCLNDBN\tCAF\tCOMMON\tMulti_base\tMulti_SNP\n')
        dGWAS = {}

        for sRow in GWAS:
            lCol = sRow.replace('\n', '').split('\t')
            sRSID_GWAS = lCol[4]
            dGWAS[sRSID_GWAS] = lCol

        for sRow in SNP:                      # Clinvar 1base.
            if sRow[0] == '#': continue
            lCol = sRow.replace('\n','').split('\t')
            sChr = lCol[0]
            if sChr == 'chrMT': continue

            sRSID = lCol[7]

            try:
                print(lCol + dGWAS[sRSID])
                GWAS_anno_output.write('\t'.join(lCol + dGWAS[sRSID]) + '\n')
            except KeyError:
                GWAS_anno_output.write('\t'.join(lCol + [""]*23) + '\n')
#def END: Check_overlap_GWAS


def Check_overlap_SNP_v2(list_sParameters):

    '''
    base editing target
    4,5,6,7base in NGG
    overlap snp.
    :param dTarget_interval:
    :return:
    '''


    sOutDir_bychr    = '%s/results_bychr' % sBE_DIR
    os.makedirs(sOutDir_bychr, exist_ok=True)

    sChrID, sStrand  = list_sParameters

    dict_sStandCheck = {'+':'F', '-':'R' }

    sInFile          = '%s/dChr_Interval_NGG_CDS_%s_%s.data' % (sTEMP_DIR, sStrand, sChrID)
    InFile           = open(sInFile, 'rb')
    dict_sTargets    = pickle.load(InFile)
    InFile.close()
    print('Interval dict loaded %s %s' % (sInFile, len(dict_sTargets)))
    sInFile          = '%s/list_cCosmic_%s_%s.data'          % (sTEMP_DIR, sChrID, sANALYSIS)
    InFile           = open(sInFile, 'rb')
    list_cCos        = pickle.load(InFile)
    InFile.close()
    print('Parsed cosmic class loaded %s %s' % (sInFile, len(list_cCos)))

    sOutFile         = '%s/Results_%s_%s_%s.txt'             % (sOutDir_bychr, sANALYSIS, sStrand, sChrID)
    dict_sOutput     = {}
    for cCos in list_cCos:

        nPos           = int(cCos.sPos)

        if dict_sStandCheck[cCos.sStrand] != sStrand: continue ## Match strands
        if not dict_sTargets[nPos]:                   continue ## Skip Empty interval i.e. no target sequence

        objCheck       = sorted(dict_sTargets[nPos])
        if not check_window(objCheck, cCos):          continue ## Check overlap, target position and window

        list_sNGGInfo  = objCheck[0][2]
        nSeqStartPos   = int(list_sNGGInfo[1])
        nSeqEndPos     = int(list_sNGGInfo[2])
        sNGG_seq       = list_sNGGInfo[3][10:33]

        if sStrand == 'F':
            sNGG_seq_ext = list_sNGGInfo[3][6:36]
        else:
            sNGG_seq_ext = list_sNGGInfo[3][7:37]

        fAlleleFreq    = check_allele_freq (cCos)

        sOutput        = [sNGG_seq, cCos, str(fAlleleFreq), sNGG_seq_ext]

        sPosKey           = '%s:%s-%s' % (sChrID, nSeqStartPos, nSeqEndPos)
        if sPosKey not in dict_sOutput:
            dict_sOutput[sPosKey] = []

        dict_sOutput[sPosKey].append(sOutput)
    #loop END: cCos

    output_results (dict_sOutput, sOutFile)
#def END: Check_overlap_SNP_v2


def Check_overlap_SNP_v2_extended(sChrIDKey):

    '''
    base editing target
    4,5,6,7base in NGG
    overlap snp.
    :param dTarget_interval:
    :return:
    '''

    sExtendFile  = '%s/%s_extend.txt' % (sBE_DIR, sANALYSIS)
    list_sExtend = load_extendfile (sExtendFile)

    dict_sRun    = {}

    for list_sParameters in list_sExtend:
        sChrID, sGeneCheck, sAltTypeCheck  = list_sParameters

        if sChrID not in dict_sRun:
            dict_sRun[sChrID] = []
        dict_sRun[sChrID].append([sGeneCheck, sAltTypeCheck])
    #loop END: list_sParameters

    sStrand       = sChrIDKey[-1]
    sChrID        = sChrIDKey[:-1]

    dict_sStandCheck = {'+':'F', '-':'R' }

    sInFile          = '%s/dChr_Interval_NGG_CDS_%s_%s.data' % (sTEMP_DIR, dict_sStandCheck[sStrand], sChrID)
    InFile           = open(sInFile, 'rb')
    dict_sTargets    =  pickle.load(InFile)
    InFile.close()
    print('Interval dict loaded %s %s' % (sInFile, len(dict_sTargets)))
    sInFile          = '%s/list_cCosmic_%s_%s.data'          % (sTEMP_DIR, sChrID, sANALYSIS)
    InFile           = open(sInFile, 'rb')
    list_cCos        = pickle.load(InFile)
    InFile.close()
    print('Parsed cosmic class loaded %s %s' % (sInFile, len(list_cCos)))

    sOutDir_bychr = '%s/results_bychr_ext' % sBE_DIR
    os.makedirs(sOutDir_bychr, exist_ok=True)

    sOutFile     = '%s/Results_extend_%s_%s.txt' % (sOutDir_bychr, sANALYSIS, sChrIDKey)
    OutFile      = open(sOutFile, 'w')

    for sCheck in dict_sRun[sChrIDKey]:

        sGeneCheck, sAltTypeCheck = sCheck


        #if sGeneCheck != 'CCDC88C': continue

        sAltTypeCheckList = sAltTypeCheck.split(',')
        #print('Check', sChrID, sStrand, sGeneCheck, sAltTypeCheckList)

        for cCos in list_cCos:

            nPos               = int(cCos.sPos)

            if cCos.sStrand   != sStrand:                  continue ## Match strands
            if cCos.sGeneName.split('_')[0] != sGeneCheck: continue
            if cCos.sAltType not in sAltTypeCheckList:     continue

            #print('Found', cCos.sGeneName.split('_')[0],  cCos.sStrand, cCos.sAltType)

            if not dict_sTargets[nPos]:                   continue ## Skip Empty interval i.e. no target sequence

            objCheck       = sorted(dict_sTargets[nPos])
            if not check_window(objCheck, cCos):          continue ## Check overlap, target position and window


            list_sNGGInfo  = objCheck[0][2]
            nSeqStartPos   = int(list_sNGGInfo[1])
            nSeqEndPos     = int(list_sNGGInfo[2])

            if sStrand == '+':
                sNGG_seq       = list_sNGGInfo[3][10:33]
                sNGG_seq_ext   = list_sNGGInfo[3][6:36]

            else:
                sNGG_seq     = reverse_complement(list_sNGGInfo[3][10:33])
                sNGG_seq_ext = reverse_complement(list_sNGGInfo[3][7:37])

            sOut           = '%s%s\t%s\t%s\t%s\t%s\n' %  (sChrID, sStrand, sGeneCheck, ','.join(sAltTypeCheckList), sNGG_seq, sNGG_seq_ext)
            OutFile.write(sOut)
        #loop END: cCos
    #loop END: sCheck
#def END: Check_overlap_SNP_v2_extended


def separate_bychr_forID_update (sInFile):

    dict_nChrID = {}
    InFile      = open(sInFile, 'r', encoding='utf-8', errors='ignore')

    for sReadLine in InFile:

        if sReadLine.startswith('Gene'):continue

        list_sColumn         = sReadLine.strip('\n').split('\t')

        cCos                 = cCOSMIC()
        cCos.sCOSID          = list_sColumn[16]
        cCos.sOldMutaID      = list_sColumn[17]
        cCos.sMutaID         = list_sColumn[18]
        cCos.sAltType        = list_sColumn[19]
        cCos.sAAType         = list_sColumn[20]
        cCos.sGenicPos       = list_sColumn[25]
        if not list_sColumn[25]: continue # Skip those w/o position information
        cCos.nChrID          = list_sColumn[25].split(':')[0]
        cCos.sStrand         = list_sColumn[26]
        cCos.bSNP            = True if list_sColumn[27] == 'y' else False
        cCos.sMutaStrand     = list_sColumn[28]
        cCos.sFATHMMPred     = list_sColumn[29] if list_sColumn[26] else 'NA'
        cCos.sFATHMMScore    = list_sColumn[30] if list_sColumn[27] else 'NA'
        cCos.sMSStatus       = list_sColumn[31] if list_sColumn[27] else 'NA'
        cCos.sPMID           = list_sColumn[32]
        cCos.sIDStudy        = list_sColumn[33]
        cCos.sSampleType     = list_sColumn[34]
        cCos.sTumorOrigin    = list_sColumn[35]
        cCos.sAge            = list_sColumn[36]
        cCos.sHGVSP          = list_sColumn[37]
        cCos.HGVSC           = list_sColumn[38]
        cCos.HGVSG           = list_sColumn[39]
        sKey                 = cCos.sAltType[-3:]

        #### KEY CHECK ####
        if sKey != dict_sKEY[sANALYSIS]:continue
        ###################

        if cCos.nChrID not in dict_nChrID:
            dict_nChrID[cCos.nChrID] = []
        dict_nChrID[cCos.nChrID].append(cCos)
    #loop END: sReadLine

    for nChrID in dict_nChrID:

        sOutFile     = '%s/COSMIC_updated_dict.%s.txt' % (sTEMP_DIR, nChrID)
        OutFile      = open(sOutFile, 'w')
        dict_sOutput = {}

        for cCos in dict_nChrID[nChrID]:

            if cCos.sOldMutaID not in dict_sOutput:
                dict_sOutput[cCos.sOldMutaID] = []
            dict_sOutput[cCos.sOldMutaID].append([cCos.sCOSID, cCos.sMutaID])
        #loop END: cCos

        for sOldID in dict_sOutput:

            if len(dict_sOutput[sOldID]) > 1:
                print(sOldID, dict_sOutput[sOldID])
        #loop END: sOldID
    #loop END: nChrID
#def END: separate_bychr_forID_update



def Check_overlap_SNP_v2_anno(list_sParameters):

    '''
    base editing target
    4,5,6,7base in NGG
    overlap snp.
    :param dTarget_interval:
    :return:
    '''

    #sOutDir_bychr    = '%s/results_bychr_anno_updated' % sBE_DIR
    sOutDir_bychr    = '%s/results_bychr_anno-Test' % sBE_DIR
    os.makedirs(sOutDir_bychr, exist_ok=True)

    sChrID, sStrand  = list_sParameters

    dict_sStandCheck = {'+':'F', '-':'R' }

    sInFile          = '%s/dChr_Interval_NGG_CDS_%s_%s.data' % (sTEMP_DIR, sStrand, sChrID)
    InFile           = open(sInFile, 'rb')
    dict_sTargets    = pickle.load(InFile)
    InFile.close()
    #print('Interval dict loaded %s %s' % (sInFile, len(dict_sTargets)))
    sInFile          = '%s/list_cCosmic_%s_%s.data'  % (sTEMP_DIR, sChrID, sANALYSIS)
    InFile           = open(sInFile, 'rb')
    list_cCos        = pickle.load(InFile)
    InFile.close()
    #print('Parsed cosmic class loaded %s %s' % (sInFile, len(list_cCos)))

    sOutFile         = '%s/Results_%s_%s_%s_anno.txt'        % (sOutDir_bychr, sANALYSIS, sStrand, sChrID)
    dict_sOutput     = {}
    for cCos in list_cCos:

        nPos           = int(cCos.sPos)

        if dict_sStandCheck[cCos.sStrand] != sStrand: continue ## Match strands
        if not dict_sTargets[nPos]:                   continue ## Skip Empty interval i.e. no target sequence

        objCheck       = sorted(dict_sTargets[nPos])
        if not check_window(objCheck, cCos):          continue ## Check overlap, target position and window


        list_sNGGInfo  = objCheck[0][2]
        nSeqStartPos   = int(list_sNGGInfo[1])
        nSeqEndPos     = int(list_sNGGInfo[2])
        sNGG_seq       = list_sNGGInfo[3][10:33]

        if sStrand == 'F':
            sNGG_seq_ext = list_sNGGInfo[3][6:36]
        else:
            sNGG_seq_ext = list_sNGGInfo[3][7:37]

        fAlleleFreq    = check_allele_freq (cCos)

        sOutput        = [sNGG_seq, cCos, str(fAlleleFreq), sNGG_seq_ext]

        sPosKey           = '%s:%s-%s' % (sChrID, nSeqStartPos, nSeqEndPos)
        if sPosKey not in dict_sOutput:
            dict_sOutput[sPosKey] = []

        dict_sOutput[sPosKey].append(sOutput)
    #loop END: cCos

    output_results_annoformat (dict_sOutput, sOutFile)
#def END: Check_overlap_SNP_v2_anno


def load_extendfile(sInFile):

    InFile       = open(sInFile, 'r')
    list_sOutput = []

    for sReadLine in InFile:

        list_sColumn = sReadLine.strip('\n').split('\t')

        list_sOutput.append(list_sColumn)

    #loop END: sReadLine

    return list_sOutput
#def END: load_extendfile


def check_window(objCheck, cCos):

    list_sNGGInfo  = objCheck[0][2]
    nSeqStartPos   = int(list_sNGGInfo[1])
    nSeqEndPos     = int(list_sNGGInfo[2])

    if cCos.sStrand == '+':
        nWinStart = nSeqStartPos + 3
        nWinEnd   = nSeqStartPos + (nWIN_SIZE - 1)
    else:
        nWinStart = nSeqEndPos - (nWIN_SIZE - 1)
        nWinEnd   = nSeqEndPos - 3

    #print(cCos.sGenicPos, cCos.sStrand, cCos.sGeneName)
    #print(sNGG_seq)
    #print(nSeqStartPos, nSeqEndPos)
    #print(nWinStart, nWinEnd)

    if nWinStart <= int(cCos.sPos) <= nWinEnd: return True
    else: return False
#def END: check_window


def check_allele_freq (cCos):

    sVCFFile   = '%s/gnomAD_exome_vcf_bychr/gnomad.exomes.%s.vcf.gz' % (sINPUT_DIR, cCos.sChrID)
    sScript    = 'tabix %s %s'                                       % (sVCFFile, cCos.sGenicPos)
    sStdOut    = sp.Popen(sScript, stdout=sp.PIPE, shell=True).stdout
    list_cVCF  = parse_vcf_stdout2(sStdOut)
    #print(sScript, cCos.sAltType, cCos.sAAType, cCos.sStrand, cCos.sGeneName)

    list_sOutput = []
    for cVCF in list_cVCF:
        #print(cVCF.nPos, cVCF.sDBSNP_ID, cVCF.sRefNuc,cVCF.sAltNuc, cVCF.fAlleleFreq, cCos.sStrand)
        if sANALYSIS == 'BE':
            if cCos.sStrand == '+':
                if cVCF.sRefNuc == 'C' and cVCF.sAltNuc == 'T':
                    list_sOutput.append(cVCF.fAlleleFreq)
            else:
                if cVCF.sRefNuc == 'G' and cVCF.sAltNuc == 'A':
                    list_sOutput.append(cVCF.fAlleleFreq)
        else:
            if cCos.sStrand == '+':
                if cVCF.sRefNuc == 'A' and cVCF.sAltNuc == 'G':
                    list_sOutput.append(cVCF.fAlleleFreq)
            else:
                if cVCF.sRefNuc == 'T' and cVCF.sAltNuc == 'C':
                    list_sOutput.append(cVCF.fAlleleFreq)
    #loop END: cVCF

    if not list_sOutput: return 'NA'
    else:                return np.mean(list_sOutput)
#def END: check_allele_freqq


def output_results (dict_sOutput, sOutFile):

    OutFile = open(sOutFile, 'w')
    sHeader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
              % ('ChrID', 'GeneSym', 'StartPos_Seq','EndPos_Seq','TargetSequence','MutationPos','gnomAD_AF','AltType','AltInfo',
                 'PrimarySite','SubSite','PrimaryHist','SubHist', 'FATHMM')
    OutFile.write(sHeader)

    for sPosKey in dict_sOutput:  #Key = chr21:43962547-43962569  #Value = [sTargetSeq, cCos, fAlleleFreq]

        sChrID, sPosition  = sPosKey.split(':')
        nStartPos, nEndPos = [int(sPos) for sPos in sPosition.split('-')]

        if len(dict_sOutput[sPosKey]) == 1:     ## Single mutation in target sequence

            cCos               = dict_sOutput[sPosKey][0][1]
            fAF                = dict_sOutput[sPosKey][0][2]
            nPos               = int(cCos.sPos)
            sGeneSym           = cCos.sGeneName.split('_')[0]
            sStrand            = cCos.sStrand
            sAltType           = cCos.sAltType
            sAltInfo           = cCos.sMutaDescri.replace(' ','')
            sPriSite           = cCos.sPriSite
            sSubSite           = cCos.sSiteSub1
            sPriHist           = cCos.sPriHist
            sSubHist           = cCos.sHistSub1
            sFATHMM            = cCos.sFATHMMPred
            sFormatSeq         = format_seq(nStartPos, nEndPos, [nPos], dict_sOutput[sPosKey][0][0], sStrand)
            sFormatSeq_ext     = dict_sOutput[sPosKey][0][3]

            sOut =  '%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    % (sChrID, sStrand, sGeneSym, nStartPos, nEndPos, sFormatSeq, sFormatSeq_ext, nPos, fAF,
                       sAltType, sAltInfo, sPriSite, sSubSite, sPriHist, sSubHist, sFATHMM)
            OutFile.write(sOut)
        else:                                   ## Multi mutation in target sequence

            list_cCos       = []
            list_sSeq       = []
            list_sSeq_ext   = []
            list_fAF        = []
            list_nPos       = []
            list_sAltType   = []
            list_sAltInfo   = []
            list_sPriSite   = []
            list_sSubSite   = []
            list_sPriHist   = []
            list_sSubHist   = []
            list_sFATHMM    = []

            for list_sInfo in dict_sOutput[sPosKey]:
                sTargetSeq, cCos, fAlleleFreq, sTargetSeq_ext = list_sInfo
                list_cCos.append(cCos)
                list_fAF.append(fAlleleFreq)
                list_sSeq.append(sTargetSeq)
                list_sSeq_ext.append(sTargetSeq_ext)
                list_nPos.append(int(cCos.sPos))
                list_sAltType.append(cCos.sAltType)
                list_sAltInfo.append(cCos.sMutaDescri.replace(' ',''))
                list_sPriSite.append(cCos.sPriSite)
                list_sSubSite.append(cCos.sSiteSub1)
                list_sPriHist.append(cCos.sPriHist)
                list_sSubHist.append(cCos.sHistSub1)
                list_sFATHMM.append(cCos.sFATHMMPred)
            #loop END: list_sInfo

            cRepCos    = list_cCos[0]
            sGeneSym   = cRepCos.sGeneName.split('_')[0]
            sStrand    = cRepCos.sStrand
            sTargetSeq     = list(set(list_sSeq))[0]
            sTargetSeq_ext = list(set(list_sSeq_ext))[0]
            sTargetSeq_ext = sTargetSeq_ext

            list_nPos  = list(set(list_nPos))

            sFormatSeq = format_seq(nStartPos, nEndPos, list_nPos, sTargetSeq, cRepCos.sStrand)

            list_fAF   = [float(fFreq) for fFreq in list_fAF if fFreq != 'NA']
            fAF        = 'NA' if not list_fAF else sum(list_fAF) / len(list_fAF)

            nPos       = '|'.join([str(nPos) for nPos in list_nPos])
            sAltType   = ','.join(sorted(list(set(list_sAltType)))[:len(list_nPos)])
            sAltInfo   = ','.join(list(set(list_sAltInfo)))
            sPriSite   = ','.join(list(set(list_sPriSite)))
            sSubSite   = ','.join(list(set(list_sSubSite)))
            sPriHist   = ','.join(list(set(list_sPriHist)))
            sSubHist   = ','.join(list(set(list_sSubHist)))
            sFATHMM    = ','.join(list(set(list_sFATHMM)))

            sOut =  '%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                        % (sChrID, sStrand, sGeneSym, nStartPos, nEndPos, sFormatSeq, sTargetSeq_ext, nPos, fAF,
                          sAltType, sAltInfo, sPriSite, sSubSite, sPriHist, sSubHist, sFATHMM)

            OutFile.write(sOut)
        #if END:
    #loop END: sPosKey
    OutFile.close()
#def END: output_results


def output_results_annoformat (dict_sOutput, sOutFile):

    OutFile = open(sOutFile, 'w')

    sHeader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
              % ('ChrID', 'StartPos_mut','EndPos_mut', 'sRefNuc', 'sAltNuc', 'Strand',
                 'GeneSym', 'SeqStart', 'SeqEnd', 'sFormatSeq', 'sFormatSeq_ext', 'gnomadAF',
                 'AltType', 'AAType', 'AltInfo', 'PrimarySite', 'SubSite', 'PrimaryHisto', 'SubHisto',
                 'MutationID', 'HGNCID', 'FATHMMPred', 'FATHMMScore', 'SomaticMutaStatus')

    OutFile.write(sHeader)

    for sPosKey in dict_sOutput:  #Key = chr21:43962547-43962569  #Value = [sTargetSeq, cCos, fAlleleFreq, sTargetSeq_ext]

        sChrID, sPosition   = sPosKey.split(':')
        nChrID              = sChrID.replace('chr','')
        nSeqStarts, nSeqEnd = [int(sPos) for sPos in sPosition.split('-')]

        if len(dict_sOutput[sPosKey]) == 1:     ## Single mutation in target sequence

            cCos               = dict_sOutput[sPosKey][0][1]
            fAF                = dict_sOutput[sPosKey][0][2]
            nPos               = int(cCos.sPos)
            sGeneSym           = cCos.sGeneName.split('_')[0]

            if sGeneSym != 'TP53': continue

            sStrand            = cCos.sStrand
            sRef, sAlt         = dict_sALTKEYS[sANALYSIS][sStrand]
            sAltType           = cCos.sAltType
            sAAType            = cCos.sAAType
            sAltInfo           = cCos.sMutaDescri.replace(' ','')
            sPriSite           = cCos.sPriSite
            sSubSite           = cCos.sSiteSub1
            sPriHist           = cCos.sPriHist
            sSubHist           = cCos.sHistSub1

            sMutaID            = cCos.sMutaID
            sHGNCID            = cCos.sHGNCID
            sFATHMMPred        = cCos.sFATHMMPred
            sFATHMMScore       = cCos.sFATHMMScore
            sMSStatus          = cCos.sMSStatus

            sTargetSeq         = dict_sOutput[sPosKey][0][0]
            sFormatSeq         = format_seq(nSeqStarts, nSeqEnd, [nPos], sTargetSeq, sStrand)

            if sStrand == '+':
                sTargetSeq_ext     = dict_sOutput[sPosKey][0][3]
                sFormatSeq_ext     = sTargetSeq_ext[:4] + sFormatSeq + sTargetSeq_ext[-3:]
            else:
                sTargetSeq_ext = reverse_complement(dict_sOutput[sPosKey][0][3])
                sFormatSeq_ext = sTargetSeq_ext[:4] + sFormatSeq + sTargetSeq_ext[-3:]

            sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                      % (nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd,
                       sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAAType, sAltInfo, sPriSite, sSubSite,
                       sPriHist, sSubHist, sMutaID, sHGNCID, sFATHMMPred, sFATHMMScore, sMSStatus)

            OutFile.write(sOut)
        else:                                   ## Multi mutation in target sequence

            list_cCos       = []
            list_sSeq       = []
            list_sSeq_ext   = []
            list_fAF        = []
            list_nPos       = []
            list_sAltType   = []
            list_sAAType    = []
            list_sAltInfo   = []
            list_sPriSite   = []
            list_sSubSite   = []
            list_sPriHist   = []
            list_sSubHist   = []

            list_sMutaID      = []
            list_sHGNCID      = []
            list_sFATHMMPred  = []
            list_sFATHMMScore = []
            list_sMSStatus    = []

            for list_sInfo in dict_sOutput[sPosKey]:
                sTargetSeq, cCos, fAlleleFreq, sTargetSeq_ext = list_sInfo
                list_cCos.append(cCos)
                list_fAF.append(fAlleleFreq)
                list_sSeq.append(sTargetSeq)
                list_sSeq_ext.append(sTargetSeq_ext)
                list_nPos.append(int(cCos.sPos))
                list_sAltType.append(cCos.sAltType)
                list_sAAType.append(cCos.sAAType)
                list_sAltInfo.append(cCos.sMutaDescri.replace(' ',''))
                list_sPriSite.append(cCos.sPriSite)
                list_sSubSite.append(cCos.sSiteSub1)
                list_sPriHist.append(cCos.sPriHist)
                list_sSubHist.append(cCos.sHistSub1)
                list_sMutaID.append(cCos.sMutaID)
                list_sHGNCID.append(cCos.sHGNCID)
                list_sFATHMMPred.append(cCos.sFATHMMPred)
                list_sFATHMMScore.append(cCos.sFATHMMScore)
                list_sMSStatus.append(cCos.sMSStatus)
            #loop END: list_sInfo

            cRepCos        = list_cCos[0]
            sGeneSym       = cRepCos.sGeneName.split('_')[0]

            if sGeneSym != 'TP53': continue


            sStrand        = cRepCos.sStrand
            sRef, sAlt     = dict_sALTKEYS[sANALYSIS][sStrand]
            sTargetSeq     = list(set(list_sSeq))[0]
            sFormatSeq     = format_seq(nSeqStarts, nSeqEnd, list_nPos, sTargetSeq, cRepCos.sStrand)

            if sStrand == '+':
                sTargetSeq_ext = list(set(list_sSeq_ext))[0]
                sFormatSeq_ext = sTargetSeq_ext[:4] + sFormatSeq + sTargetSeq_ext[-3:]
            else:
                sTargetSeq_ext = reverse_complement(list(set(list_sSeq_ext))[0])
                sFormatSeq_ext = sTargetSeq_ext[:4] + sFormatSeq + sTargetSeq_ext[-3:]


            list_nPos  = list(set(list_nPos))

            list_fAF   = [float(fFreq) for fFreq in list_fAF if fFreq != 'NA']
            fAF        = 'NA' if not list_fAF else sum(list_fAF) / len(list_fAF)

            nPos       = '|'.join([str(nPos) for nPos in list_nPos])
            sAltType   = ','.join(sorted(list(set(list_sAltType)))[:len(list_nPos)])
            sAAType    = ','.join(sorted(list(set(list_sAAType)))[:len(list_nPos)])
            sAltInfo   = ','.join(list(set(list_sAltInfo)))
            sPriSite   = ','.join(list(set(list_sPriSite)))
            sSubSite   = ','.join(list(set(list_sSubSite)))
            sPriHist   = ','.join(list(set(list_sPriHist)))
            sSubHist   = ','.join(list(set(list_sSubHist)))

            sMutaID      = ','.join(list(set(list_sMutaID)))
            sHGNCID      = ','.join(list(set(list_sHGNCID)))
            sFATHMMPred  = ','.join(list(set(list_sFATHMMPred)))
            sFATHMMScore = ','.join(list(set(list_sFATHMMScore)))
            sMSStatus    = ','.join(list(set(list_sMSStatus)))

            sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                      % (nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd,
                       sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAAType, sAltInfo, sPriSite, sSubSite,
                       sPriHist, sSubHist, sMutaID, sHGNCID, sFATHMMPred, sFATHMMScore, sMSStatus)


            OutFile.write(sOut)
        #if END:
    #loop END: sPosKey
    OutFile.close()
#def END: output_results_annoformat


def format_seq (nStartPos, nEndPos, list_nAltPos, sTargetSeq, sStrand):

    list_nAltIndex = []
    for nAltPos in list_nAltPos:
        if sStrand == '+':
            nWinStart  = nStartPos + 3
            nWinEnd    = nStartPos + nWIN_SIZE

            nWinIndexS = nWinStart - nStartPos
            nWinIndexE = nWinEnd   - nStartPos
            nAltIndex  = nAltPos - nStartPos

        else:
            nWinStart = nEndPos - nWIN_SIZE
            nWinEnd   = nEndPos - 3

            nWinIndexS = nWinStart - nStartPos + 1
            nWinIndexE = nWinEnd   - nStartPos + 1
            nAltIndex  = nAltPos   - nStartPos
        #if END: sStrand

        list_nAltIndex.append(nAltIndex)

    #loop END: list_nAltIndex

    #print(nStartPos, nEndPos,nAltPos)
    #print(nWinStart, nWinEnd)
    #print(nWinIndexS, nWinIndexE, nAltIndex)

    list_sFormatSeq = []
    for i in range(len(sTargetSeq)):
        #print(i, sTargetSeq[i])
        if i in list_nAltIndex:  list_sFormatSeq.append(sTargetSeq[i].lower())
        else: list_sFormatSeq.append(sTargetSeq[i])
    #loop END: i

    sFormatSeq = ''.join([sStr for sStr in list_sFormatSeq])
    if sStrand == '+': return sFormatSeq
    else:              return reverse_complement(sFormatSeq)
    return sFormatSeq
#def END: format_seq


def  combine_sort_output (list_sParameters, sFileTag):

    sOutDir_bychr  = '%s/results_bychr_anno' % sBE_DIR
    sOutDir        = '%s/final'              % sBE_DIR
    sLogDir        = '%s/log'                % sBE_DIR
    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    list_sFiles    = []
    for sParameters in list_sParameters:
        sChrID, sStrand = sParameters
        sInFile         = '%s/Results_%s_%s_%s_anno.txt' % (sOutDir_bychr, sANALYSIS, sStrand, sChrID)
        list_sFiles.append(sInFile)
    #loop END: sParameters

    sOutFile  = '%s/%s.temp.txt'   % (sOutDir, sFileTag)
    sOutFile2 = '%s/%s.sorted.txt' % (sOutDir, sFileTag)
    sOutFile3 = '%s/%s.final.txt'  % (sOutDir, sFileTag)

    sScript   = 'cat %s > %s;'                               % (' '.join(list_sFiles), sOutFile)
    sScript  += 'sort -k1,1g -k3,3g  -S 12G -T %s %s -o %s;' % (sLogDir, sOutFile, sOutFile2)
    sScript  += 'grep -v ChrID %s > %s; '                    % (sOutFile2, sOutFile3)
    sScript  += 'rm %s %s'                                   % (sOutFile, sOutFile2)

    os.system(sScript)
#def END: combine_sort_output


def basic_stats_final_output (sFileTag):

    sInDir     = '%s/final'         % sBE_DIR
    sInFile    = '%s/%s.final.txt'  % (sInDir, sFileTag)
    InFile     = open(sInFile, 'r')
    dict_nCnts = {}
    nCnt       = 0
    for sReadLine in InFile:

        list_sColumn = sReadLine.strip('\n').split('\t')
        sAltType     = list_sColumn[14]
        nCnt        += 1

        if sAltType not in dict_nCnts:
            dict_nCnts[sAltType] = 0
        dict_nCnts[sAltType] += 1

    #loop END: sReadLine
    InFile.close()

    print('Total', nCnt)
    for sAltType in dict_nCnts:
        print(sAltType, dict_nCnts[sAltType])
#def END: basic_stats_final_output


def update_COSMIC_info (sFileTag):
    sOutDir       = '%s/final'                 % sBE_DIR
    sInFile       = '%s/%s.final.txt'          % (sOutDir, sFileTag)
    sOutFile      = '%s/%s.final_Updated.txt'  % (sOutDir, sFileTag)
    OutFile       = open(sOutFile, 'w')

    dict_sOutput  = load_final_file (sInFile)
    dict_sNewCnt  = {'UNCHANGED':0, 'UPDATED':0}
    for nChrID in dict_sOutput:
        list_sData       = dict_sOutput[nChrID]
        sInFile          = '%s/list_cCosmic_updated_chr%s_%s.data'  % (sTEMP_DIR, nChrID, sANALYSIS)
        InFile           = open(sInFile, 'rb')
        list_cCos        = pickle.load(InFile)
        InFile.close()
        print(nChrID, len(list_sData), len(list_cCos))

        dict_cCos = {}
        for cCos in list_cCos:

            sKey = cCos.sOldMutaID

            #if sKey in ['COSM984916','COSM984914','COSM10777','COSM1717141','COSM3701291','COSM984917']:
                #print(cCos.sGeneName, cCos.sOldMutaID, cCos.sCOSID, cCos.sAAType, cCos.sAltType)

            if sKey not in dict_cCos:
                dict_cCos[sKey] = []
            dict_cCos[sKey].append(cCos)
        #loop END: cCos

        for sData in list_sData:
            nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd, \
            sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAAType, sAltInfo, sPriSite, \
            sSubSite, sPriHist, sSubHist, sMutaID, sHGNCID, sFATHMMPred, sFATHMMScore, sMSStatus = sData

            if ',' in sMutaID: list_sMutaID = sMutaID.split(',')
            else: list_sMutaID = [sMutaID]

            list_sCOSID   = []
            list_sAAType  = []
            list_sAltType = []

            for sMutaID in list_sMutaID:

                try: list_cCos = dict_cCos[sMutaID]
                except KeyError: continue

                for cCos in list_cCos:

                    if cCos.sGeneName != sGeneSym: continue
                    list_sCOSID.append(cCos.sCOSID)
                    list_sAAType.append(cCos.sAAType)
                    list_sAltType.append(cCos.sAltType)
                #loop END: sMutaID
            #loop END: sMutaID

            sNewID      = ','.join(list(set(list_sCOSID)))
            sNewAAType  = ','.join(list(set(list_sAAType)))
            sNewAltType = ','.join(list(set(list_sAltType)))

            if sNewID and sNewAAType and sNewAltType:
                sOut = '\t'.join([nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd, \
                                  sFormatSeq, sFormatSeq_ext, fAF, sNewAltType, sNewAAType, sAltInfo, sPriSite, \
                                sSubSite, sPriHist, sSubHist, sNewID, sHGNCID, sFATHMMPred, sFATHMMScore, sMSStatus])
                #print('NEW', sOut)
                dict_sNewCnt['UPDATED'] += 1

            else:
                sOut = '\t'.join(sData)
                #print('UNCHANGED', sOut)
                dict_sNewCnt['UNCHANGED'] += 1
            #if END:
            OutFile.write('%s\n' % sOut)
        #loop END: sData
    #loop END: nChrID
    OutFile.close()
    for sNew in dict_sNewCnt:
        print(sNew, dict_sNewCnt[sNew])
#def END: update_COSMIC_info


def load_final_file (sInFile):

    dict_sOutput = {}
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:

        list_sColumn = sReadLine.strip('\n').split('\t')
        nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd, \
        sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAAType, sAltInfo, sPriSite, \
        sSubSite, sPriHist, sSubHist, sMutaID, sHGNCID, sFATHMMPred, sFATHMMScore, sMSStatus = list_sColumn

        if nChrID not in dict_sOutput:
            dict_sOutput[nChrID] = []
        dict_sOutput[nChrID].append(list_sColumn)

    #loop END: sReadLine
    InFile.close()

    return dict_sOutput
#def END: load_final_file


def check_clinvar_overlap (sFileTag):

    sOutDir       = '%s/final'        % sBE_DIR
    sInFile       = '%s/%s.final.txt' % (sOutDir, sFileTag)
    sClinvarFile  = '%s/clinvar/clinvar_20210328.vcf.gz' % sREF_DIR

    sOutFile      = '%s/%s.final_ClinVarMatch.txt'  % (sOutDir, sFileTag)
    OutFile       = open(sOutFile, 'w')
    dict_sOutput  = load_final_file (sInFile)

    for nChrID in dict_sOutput:

        list_sData = dict_sOutput[nChrID]

        for sData in list_sData:
            nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd, \
            sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAAType, sAltInfo, sPriSite, \
            sSubSite, sPriHist, sSubHist, sMutaID, sHGNCID, sFATHMMPred, sFATHMMScore, sMSStatus = sData

            if len(nPos.split('|')) > 1:
                nPos1 = nPos.split('|')[0]
                nPos2 = nPos.split('|')[1]
            else:
                nPos1 = nPos
                nPos2 = nPos

            sCmd      = 'tabix %s %s:%s-%s' % (sClinvarFile, nChrID, nPos1, nPos2)
            sStdOut   = subprocess.Popen(sCmd, stdout=subprocess.PIPE, shell=True).stdout
            list_cVCF = cVCF_parse_vcf_files_clinvar_stdout(sStdOut, sRef, sAlt)

            if list_cVCF:

                for cVCF in list_cVCF:
                    sClinvarOut = '%s,%s,%s,%s,%s,%s,%s' % (cVCF.sChrID, cVCF.nPos, cVCF.sRefNuc, cVCF.sAltNuc, cVCF.fQual, cVCF.sFilter, cVCF.sInfo )
                    sOut        = '\t'.join([nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd, \
                                      sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAAType, sAltInfo, sPriSite, \
                                      sSubSite, sPriHist, sSubHist, sMutaID, sHGNCID, sFATHMMPred, sFATHMMScore,
                                      sMSStatus, sClinvarOut])
                    OutFile.write('%s\n' % sOut)

            else:
                sOut = '\t'.join(sData)
                OutFile.write('%s\n' % sOut)
            # if END:

        #loop END: sData
    #loop END: nChrID
    OutFile.close()
#def END: check_clinvar_overlap


def check_clinvar_overlap_v2 (sFileTag):

    sInFile       = '%s/20210406_DoubleCheckClinvar%s.txt' % (sINPUT_DIR, sANALYSIS)
    dict_sInput   = load_overlap_data (sInFile)

    sClinvarFile  = '%s/clinvar/clinvar_20210328.vcf.gz' % sREF_DIR

    ## COSMIC output file ##
    #sFinalDir     = '%s/final'        % sBE_DIR
    #sInFile       = '%s/%s.final.txt' % (sFinalDir, sFileTag)
    #dict_sOutput  = load_final_file (sInFile)

    sOutFile      = '%s/20210408_DoubleCheck.output%s.txt'   % (sOUTPUT_DIR, sANALYSIS)
    OutFile       = open(sOutFile, 'w')
    for sChrID in dict_sInput:

        list_sData = dict_sInput[sChrID]

        for sData in list_sData:
            sLibrary, sSG_ID, sMainID, sGuideSeq, sTarSeq, sPosKey, sGeneSym, sAA = sData

            sChrID, sPos = sPosKey.split(':')

            #if ',' not in sSG_ID: continue

            sCmd      = 'tabix %s %s:%s-%s' % (sClinvarFile, sChrID, sPos, sPos)

            sStdOut   = subprocess.Popen(sCmd, stdout=subprocess.PIPE, shell=True).stdout
            list_cVCF = cVCF_parse_vcf_files_clinvar_stdout(sStdOut)


            if list_cVCF:
                for cVCF in list_cVCF:
                    sClinvarOut = '%s,%s,%s,%s,%s,%s,%s' % (cVCF.sChrID, cVCF.nPos, cVCF.sRefNuc, cVCF.sAltNuc, cVCF.fQual, cVCF.sFilter, cVCF.sInfo )
                    sOut        = '\t'.join([sLibrary, sSG_ID, sMainID, sGuideSeq, sTarSeq, sPosKey, sGeneSym, sAA, sClinvarOut])

                    #print(sOut)
                    OutFile.write('%s\n' % sOut)

            #else:
                #sOut = '\t'.join(sData)
                #print(sOut)
                #OutFile.write('%s\n' % sOut)
            # if END:
        #loop END: sData
    #loop END: nChrID
    OutFile.close()
#def END: check_clinvar_overlap_v2


def load_overlap_data (sInFile):

    InFile       = open(sInFile, 'r')
    dict_sOutput = {}
    for sReadLine in InFile:

        # File Format
        # sg.ID	surrogate_ID		Guide Sequence	GX19	Final_classification	Gene	mutation_pattern	Clinvar
        if sReadLine.startswith('#Library'): continue  # SKIP Information Headers
        list_sColumn = sReadLine.strip('\n').split('\t')

        sLibrary    = list_sColumn[0]
        sSG_ID      = list_sColumn[1]
        sMainID     = list_sColumn[2].split('_')[1]
        sGuideSeq   = list_sColumn[3]
        sTarSeq     = list_sColumn[4]
        sPosKey     = list_sColumn[5]
        sGeneSym    = list_sColumn[6]
        sMainCheck  = list_sColumn[7]

        #print(sLibrary, sSG_ID, sMainID, sGuideSeq, sTarSeq, sPosKey, sGeneSym, sMainCheck)

        sChrID    = sPosKey.split(':')[0]
        list_sPos = sPosKey.split(':')[1].split('|')
        list_sAA  = sSG_ID.split('_')[1].split(',')

        if sChrID not in dict_sOutput:
            dict_sOutput[sChrID] = []

        for sAA, sPos in zip(list_sAA, list_sPos):

            sNewPosKey = '%s:%s' % (sChrID, sPos)

            dict_sOutput[sChrID].append([sLibrary, sSG_ID, sMainID, sGuideSeq, sTarSeq, sNewPosKey, sGeneSym, sAA])
        #loop END: sPos
    #loop END: sReadLine

    return dict_sOutput

#def END: load_overlap_data



def cVCF_parse_vcf_files_clinvar_stdout (sStdOut):
    sRef, sAlt = dict_sKEY[sANALYSIS].split('>')
    list_sOutput = []
    for sReadLine in sStdOut:
        sReadLine = str(sReadLine, 'UTF-8').strip('\n')
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | 1       | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat

        if sReadLine.startswith('#'): continue  # SKIP Information Headers
        list_sColumn = sReadLine.strip('\n').split('\t')

        cVCF         = cVCFData()
        cVCF.sChrID  = 'chr%s' % list_sColumn[0]

        if list_sColumn[0].startswith('MT'): continue
        if list_sColumn[0].startswith('NW'): continue

        try:               cVCF.nPos = int(list_sColumn[1])
        except ValueError: continue

        cVCF.sDBSNP_ID = list_sColumn[2]
        cVCF.sRefNuc = list_sColumn[3]
        cVCF.sAltNuc = list_sColumn[4]
        cVCF.fQual = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
        cVCF.sFilter = list_sColumn[6]
        cVCF.sInfo = list_sColumn[7]

        sKey      = '%s>%s' % (cVCF.sRefNuc, cVCF.sAltNuc)
        sKeyCheck = '%s>%s' % (sRef, sAlt)
        if sKey != sKeyCheck: continue

        dict_sInfo = dict([sInfo.split('=') for sInfo in cVCF.sInfo.split(';') if len(sInfo.split('=')) == 2])

        # try: cVCF.sAlleleFreq    = float(dict_sInfo['AF_raw'])
        # except ValueError: cVCF.sAlleleFreq = np.mean([float(f) for f in  dict_sInfo['AF_raw'].split(',')])

        list_sOutput.append(cVCF)
    #loop END: sReadLine

    return list_sOutput
#def END: cVCF_parse_vcf_files_clinvar


def convert_to_vcf (sFileTag):

    sWorkDir = '%s/final'          % sBE_DIR
    sInFile  = '%s/%s.final.txt'   % (sWorkDir, sFileTag)
    InFile   = open(sInFile, 'r')

    sOutFile = '%s/%s.final.vcf'   % (sWorkDir, sFileTag)
    OutFile  = open(sOutFile, 'w')
    for sReadLine in InFile:

        #VCF Format
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

        list_sColumn   = sReadLine.replace('\n','').split('\t')
        nChrID         = list_sColumn[0]
        nPos           = list_sColumn[1]
        sID            = '-'
        sRef           = list_sColumn[3]
        sAlt           = list_sColumn[4]
        sQual          = '-'
        sFilter        = 'PASS'
        sInfo          = '|'.join(list_sColumn[5:])

        sOut =  '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                % (nChrID, nPos, sID, sRef, sAlt, sQual, sFilter, sInfo)


        #sOut =  '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        #         % (nChrID, nPos, nPos, sRef, sAlt, sStrand, sGeneSym, nSeqStarts, nSeqEnd,
        #           sFormatSeq, sFormatSeq_ext, fAF, sAltType, sAltInfo, sPriSite, sSubSite,
        #           sPriHist, sSubHist, sFATHMM)
        OutFile.write(sOut)
    #loop END: sReadLine
    InFile.close()
    OutFile.close()
#def END: convert_to_vcf


def annotate_annovar (sFileTag):

    sInDir   = '%s/final'              % sBE_DIR
    sInFile  = '%s/%s.final.vcf'       % (sInDir, sFileTag)

    sOutDir  = '%s/final_anno'         % sBE_DIR
    sOutTag  = '%s/%s.final.anno'      % (sOutDir, sFileTag)

    sLogDir  = '%s/log'                % sBE_DIR
    sLogFile = '%s/log_%s.txt'         % (sLogDir, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)
    os.makedirs(sOutDir, exist_ok=True)

    # Annotate with ANNOVAR
    sScript     = '%s -filter '         % sANNOVAR
    sScript    += '-dbtype ALL.sites.2015_08 '
    sScript    += '-buildver hg19 '
    sScript    += '%s '                 % sInFile
    sScript    += '%s '                 % sDB_DIR
    sScript    += '--outfile %s '       % sOutDir
    sScript    += '-otherinfo; '

    sScript      = '%s -format vcf %s ' % (sAVINPUT, sInFile)
    sScript     += '--outfile %s.test ' % sOutTag
    sScript     += '-allsample '
    sScript     += '-withfreq; '

    #os.system(sScript)
    print(sScript)
#def END: annotate_annovar


def parse_genelist (sInFile):
    return [sReadLine.replace('\n','').upper() for sReadLine in open(sInFile, 'r')]
#def END: parse_genelist


def get_gene_info (sWorkDir, list_sTargetGenes):

    sOutDir      = '%s/OutputPerGene' % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)

    list_sKeys   = ['cellline', 'chr', 'start', 'end', 'symbol', 'sequence', 'strand', 'pubmed',
                    'cas', 'screentype', 'condition', 'effect',  'ensg', 'log2fc']

    sHeader       = '\t'.join(list_sKeys)

    for sGene in list_sTargetGenes:

        print('Processing %s' % sGene)

        sOutFile = '%s/%s.output.txt' % (sOutDir, sGene)
        OutFile  = open(sOutFile, 'w')
        OutFile.write('%s\n' % sHeader)

        sScript  = 'curl -s -H "Content-Type: application/json" '
        sScript += '-X POST -d \'{"query":"%s"}\' ' % sGene
        sScript += 'http://genomecrispr.dkfz.de/api/sgrnas/symbol; '

        sStdOut    = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout
        list_sJSON = json.load(sStdOut)


        for dict_sInfo in list_sJSON:

            list_sInfo = ['%s' % dict_sInfo[sKey] for sKey in list_sKeys]

            sOut       = '\t'.join(list_sInfo)
            OutFile.write('%s\n' % sOut)
        #loop END: dict_sInfo
        OutFile.close()

        print('Processing %s.............Done' % sGene)
    #loop END: sGene
#def END: get_gene_info


def get_lines_wTargetGenes (list_sTargetGenes, sGeneInfoFile, sOutputFile):

    if not os.path.isfile(sGeneInfoFile): sys.exit('File Not Found %s' % sGeneInfoFile)

    InFile  = open(sGeneInfoFile, 'r')
    OutFile = open(sOutputFile, 'w')

    for sReadLine in InFile:

        if sReadLine.startswith('start'): continue

        list_sColumn = sReadLine.replace('\n','').split(',')

        sGeneSym     = list_sColumn[8]

        if sGeneSym in list_sTargetGenes: OutFile.write(sReadLine)

    #loop END: sReadLine
    InFile.close()
    OutFile.close()
#def END: get_lines_wTargetGenes


def temp_lookup ():

    InFile1 = '%s/lookup/reference.txt' % sBASE_DIR
    InFile2 = '%s/lookup/position.txt'  % sBASE_DIR

    list_sRef = [sReadLine.replace('\n', '') for sReadLine in open(InFile1)]
    list_sPos = [sReadLine.replace('\n', '') for sReadLine in open(InFile2)]

    print(list_sRef[:10])
    print(list_sPos[:10])

    dict_sPos = {sPos : i for i,sPos in enumerate(list_sPos)}


    sOutFile  = '%s/lookup/ReferenceCheck.txt'  % sBASE_DIR
    OutFile   = open(sOutFile, 'w')
    for sRef in list_sRef:

        try:             sOut = '%s\t%s\n' % (sRef, dict_sPos[sRef])
        except KeyError: sOut = '%s\t%s\n' % (sRef, 'NotFound')

        OutFile.write(sOut)

    #loop END: sRef
    OutFile.close()

#def END: temp_lookup


def main():

    cores = 24
    p     = mp.Pool(processes=cores)

    list_sParameters = [[sChrID, sStrand] for sChrID in list_sCHRIDs for sStrand in ['F','R']]

    ## Analysis Settings Conditions
    print(cyan('**********Analysis Conditions**********',1))
    print(cyan('Study and WinSize\t%s-%s' % (sANALYSIS, nWIN_SIZE),1))
    print(cyan('***************************************',1))
    print ('program start: %s' % datetime.now())

    ### STEP 1 ###
    #print 'Set interval tree'
    dict_chr_interval = set_intervalTree(Path.sCCDS)
    ###############

    ### STEP  2 ###
    #print 'Filter NGG 20mer'
    p.map_async(filter_NGG_20mer, [[dict_chr_interval, 'F'],[dict_chr_interval, 'R']]).get()
    ###############

    ### STEP 3 ###
    #print ('Set interval tree base editing target')
    p.map_async(set_interval_tree_v2, list_sParameters).get()
    ###############

    ### STEP 4 ###
    #print ('gnomAD annotation')
    p.map_async(make_GNOMAD_by_chr, list_nCHRIDs).get()
    ##############

    ### STEP 5 ###
    #print ('Parse cosmic data')
    cCos_parse_cosmic_consensus('chrY')
    p.map_async(cCos_parse_cosmic_consensus, list_sCHRIDs).get()
    ##############

    ### STEP 6 ###
    #print('Check Target Sequence and Mutations')
    sNewCosmicFile  = '%s/COSMIC/cosmic_mutations_hg19_updated.tsv' % sREF_DIR
    separate_bychr_forID_update (sNewCosmicFile)
    #sys.exit()

    #Check_overlap_SNP_v2_anno(['chr17','R'])
    p.map_async(Check_overlap_SNP_v2_anno, list_sParameters).get()
    ##############kjo

    ### STEP 7 ###
    #sFileTag = '200916_wAA-COSMIC_gnomAD_results_multi_CCDS_%s' % sANALYSIS
    sFileTag = '210325_wAA-COSMIC_gnomAD_results_multi_CCDS_%s' % sANALYSIS ## TEST
    #print('Combine and Sort Output')
    #combine_sort_output (list_sParameters, sFileTag)
    #basic_stats_final_output (sFileTag)
    #update_COSMIC_info (sFileTag)
    #check_clinvar_overlap (sFileTag)
    #check_clinvar_overlap_v2 (sFileTag)
    ##############

    ### STEP 7 ###
    #print('Annotated with ANNOVAR')
    #convert_to_vcf (sFileTag)
    annotate_annovar(sFileTag)
    ##############

    '''

    ### Short Analyses ###
    sWorkDir           = '%s/Input/genomeCRISPR' % sBASE_DIR
    #get_gene_info (sWorkDir, list_sTargetGenes)

    sWorkDir           = '%s/Input/genomeCRISPR/' % sBASE_DIR
    sTargetGenesFile   = '%s/targetgenes.txt'     % sWorkDir
    sGeneInfoFile      = '%s/GenomeCRISPR_full05112017.csv' % sWorkDir
    sGeneInfoOutFile   = '%s/GenomeCRISPR_full05112017_New.output.txt' % sWorkDir

    list_sTargetGenes  = parse_genelist(sTargetGenesFile)
    print('Target Genes', len(list_sTargetGenes))

    get_lines_wTargetGenes (list_sTargetGenes, sGeneInfoFile, sGeneInfoOutFile)

    print('program end: %s' % datetime.now())

    ## Short Analyses for Myungae - Genome Sequence Inquiry
    sWorkDir        = '%s/GenomeInquiry'        % sBASE_DIR
    nBufferSize     = 0
    sLabel          = 'forHKK3'
    list_sGenome    = ['hg38','hg19']
    list_sGenome    = ['hg38']

    sOutFile        = '%s/GenomeTargetSequences_%s_Buff%s.txt' % (sWorkDir, sLabel, nBufferSize)
    OutFile         = open(sOutFile, 'w')

    for sGenome in list_sGenome:
        sGenomeFile = '/data/reference_genome/%s/%s.fa' % (sGenome, sGenome)
        sTargetFile = '%s/TargetLocal%s_buff%s_%s.txt'  % (sWorkDir,sLabel, nBufferSize, sGenome)
        list_sTarget = [sTarget.replace('\n', '').split('\t') for sTarget in open(sTargetFile, 'r')]

        for sID, sTarget in list_sTarget:
            sChrID, sTarPos = sTarget.split(':')
            sChrID = sChrID.replace('Chr', 'chr')
            if not sChrID.startswith('chr'): sChrID = 'chr%s' % sChrID

            nStart, nEnd    = [int(sPos) for sPos in sTarPos.split('-')]

            nTarStart       = nStart - nBufferSize
            nTarEnd         = nEnd   + nBufferSize

            cGenome         = cFasta(sGenomeFile)
            try: sTarSeq    = cGenome.fetch(sChrID.replace('Chr','chr'), nTarStart-1, nTarEnd).upper()
            except  AssertionError:
                print(sID, sTarget)
                continue
            #print(sID, sChrID, nStart, nEnd, sTarSeq, len(sTarSeq), nTarEnd - nTarStart)

            sOut            = '%s\t%s\t%s:%s-%s\t%s\t%s\n' \
                              % (sID, sTarget, sChrID, nTarStart, nTarEnd,  sTarSeq, sGenome)
            OutFile.write(sOut)
        #loop END: sID, sTarget
    #loop END: sGenome
    OutFile.close()
    '''

    '''

    sExtendFile  = '%s/%s_extend.txt' % (sBE_DIR, sANALYSIS)
    list_sExtend = load_extendfile (sExtendFile)
    dict_sRun    = {}

    for list_sParameters in list_sExtend:
        sChrID, sGeneCheck, sAltTypeCheck  = list_sParameters

        if sChrID not in dict_sRun:
            dict_sRun[sChrID] = []
        dict_sRun[sChrID].append([sGeneCheck, sAltTypeCheck])
    #loop END: list_sParameters
    #Check_overlap_SNP_v2_extended('chr14-')
    p.map_async(Check_overlap_SNP_v2_extended, dict_sRun).get()

    #print 'EXAC annotation'
    #Make_EXAC_dict()

    #print('Annotate GWAS')
    #Check_overlap_SNP(dNGG_with_NGG_CDS_forward_extend, 'F', dict_sDatabase, sAnalysis)
    #Check_overlap_SNP(dNGG_with_NGG_CDS_reverse_extend, 'R', dict_sDatabase, sAnalysis)
    
    #Check_overlap_GWAS('F')
    #Check_overlap_GWAS('R')

    #print('Extend_200bp_with_NGG_20mer')
    # print('Make hg19 dict')
    # dHg19 = Make_hg19_dict(
    #Extend_200bp_with_NGG_20mer(dHg19, iRandom_selection_num)
    '''
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