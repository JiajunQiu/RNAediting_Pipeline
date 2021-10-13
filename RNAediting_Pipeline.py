import glob
import os 
from scipy.stats import norm
import re
from statsmodels.sandbox.stats.multicomp import multipletests
import pysam
from optparse import OptionParser

#Commandline parsing
disc = "RNA editing identification pipeline"
#usage = "usage: %prog [options]"+'\n\nExample:\npython RNAediting_Pipeline -i ./dataset/test_pred.txt\n'
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage,description = disc)
parser.add_option("-i", action="store", type="string", dest="bam", help="A bam file from sequence reads alignment program")
parser.add_option("-g", action="store", type="string", dest="genome", help="human genome fasta file")
parser.add_option("-s", action="store", type="string", dest="snp", help="vcf file for known SNPs")

options, args = parser.parse_args()

bam_file = options.bam
genome_file = options.genome
snp_file = options.snp

def RNA_editing_identify(bam_file)
    dirr='/'.join(bam_file.split('/')[:-1])
    if dirr=='':
        dirr='.'
    tag=bam_file.split('.bam')[0]

    bam=bam_file
    
    
    tmp_out1=dirr+'/'+'test1.bam'
    tmp_out2=dirr+'/'+'test2.bam'
    tmp_out3=dirr+'/'+'test3.bam'
    tmp_out4=dirr+'/'+'test4.bam'
  

    os.system('bam trimBam '+bam+' '+dirr+'/'+'tmp.bam -L 6 -c')
    os.system('java -Xmx15000m -jar ./picard.jar SortSam  TMP_DIR='+dirr+'/'+'temp I='+dirr+'/'+'tmp.bam O='+dirr+'/'+'tempQuerySort.bam SORT_ORDER=queryname')
#    os.system('samtools sort -n  tmp.bam -o tempQuerySort.bam')
    os.unlink(dirr+'/'+'tmp.bam')
#    os.system('samtools fixmate tempQuerySort.bam  tmp.bam')
    os.system('java -Xmx15000m -jar ./picard.jar FixMateInformation I='+dirr+'/'+'tempQuerySort.bam O='+dirr+'/'+'tmp.bam')
    os.unlink(dirr+'/'+'tempQuerySort.bam')
    os.system('java -Xmx15000m -jar ./picard.jar SortSam  TMP_DIR='+dirr+'/'+'temp I='+dirr+'/'+'tmp.bam O='+tmp_out1+' SORT_ORDER=coordinate')
    os.unlink(dirr+'/'+'tmp.bam')
    os.system('java -Xmx15000m -jar ./picard.jar MarkDuplicates  TMP_DIR='+dirr+'/'+'temp I='+tmp_out1+' O='+tmp_out2+' M='+dirr+'/'+'test.txt')
    os.system('java -Xmx15000m -jar ./picard.jar AddOrReplaceReadGroups  TMP_DIR='+dirr+'/'+'temp I='+tmp_out2+' O='+tmp_out3+' RGID=SRR RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM='+dirr.split('/')[-1])

    os.system('java -Xmx15000m -jar ./picard.jar ReorderSam TMP_DIR='+dirr+'/'+'temp I='+tmp_out3+' O='+tmp_out4+' REFERENCE= '+genome_file+' ALLOW_INCOMPLETE_DICT_CONCORDANCE=True')

    os.system('samtools index '+tmp_out4)
    
    os.system('gatk3 -Xmx32g -T  BaseRecalibrator -nct 16 -R '+genome_file+' -I '+tmp_out4+' -knownSites '+snp_file+' -o '+dirr+'/recal_data.table -U ALLOW_N_CIGAR_READS')
    os.system('gatk3 -Xmx32g -T  HaplotypeCaller -R '+genome_file+' -I '+tmp_out4+' -BQSR '+dirr+'/recal_data.table -o '+dirr+'/raw_variants.vcf  --output_mode EMIT_ALL_SITES -U ALLOW_N_CIGAR_READS --min_base_quality_score 25 --min_mapping_quality_score 20 -stand_call_conf 0 -bamout '+dirr+'/tmp.bam')     


    os.unlink(tmp_out1)
    os.unlink(tmp_out2)
    os.unlink(tmp_out3)
    
    f=dirr+'/raw_variants.vcf'

    head=[]
    for l in open(f):
        l=l.rstrip()
        if '#' in l:
            head.append(l)
        else:
            break



    #filter1
    os.system('gatk3 -Xmx32g -T SelectVariants -R '+genome_file+' -V '+dirr+'/raw_variants.vcf'+' --discordance '+snp_file+' -o '+dirr+'/filter1_variants.vcf')
    

    #filter2
    out=open(dirr+'/filter2_variants.vcf','w')
    p_values=[]
    data=[]
    dp=[]
    alt=[]
    inde=[]
    inde_num={}
    for l in open(dirr+'/filter1_variants.vcf'):
        l=l.rstrip()
        if l[:2]!='ch':
            print(l,file=out)
        elif 'ReadPosRankSum' in l and len(l.split('\t')[3])==1 and len(l.split('\t')[4])==1 and l.split('\t')[-1]!='./.':
            m = re.search(r'ReadPosRankSum=(.*?);',l)
            z_score = float(m.group(1))
            p = norm.sf(z_score)*2
            p_values.append(p)
            data.append(l)
            temp=l.split('\t')[-1]
            inde.append((l.split('\t')[0],l.split('\t')[1]))
            dp.append(float(temp.split(':')[2]))
            alt.append(float(temp.split(':')[1].split(',')[1]))
            if (l.split('\t')[0],l.split('\t')[1]) in inde_num:
                inde_num[(l.split('\t')[0],l.split('\t')[1])]+=1
            else:
                inde_num[(l.split('\t')[0],l.split('\t')[1])]=1
    
    p_adjusted = multipletests(p_values, method='bonferroni')

    for n in range(len(p_values)):
        if p_values[n]>=0.01 and alt[n]/dp[n]>=0.1 and inde_num[inde[n]]==1:
            print(data[n],file=out)


    


    #filter3
    os.system('bam splitBam -i '+dirr+'/tmp.bam'+' -o '+dirr+'/splited')

    os.system('bedtools intersect -u -wa -a '+dirr+'/splited.SRR.bam'+' -b '+dirr+'/filter2_variants.vcf'+' > '+dirr+'/tmp2.bam')
    
    
   


    os.unlink(dirr+'/splited.SRR.bam')
    os.unlink(dirr+'/splited.ArtificialHaplotype.bam')
    
    



    os.system('samtools view -b -f 0x40 '+dirr+'/tmp2.bam'+' -o '+dirr+'/tmp2_1.bam')



    os.system('bedtools intersect -wb -a '+dirr+'/filter2_variants.vcf'+' -b '+dirr+'/tmp2_1.bam'+' > '+dirr+'/variant_reads')

    os.system('bedtools bamtofastq -i '+dirr+'/tmp2_1.bam -fq '+dirr+'/tmp.fa')
    
    var_dp={}
    reads_map_pos={}
    pos_map_reads={}
    count={}
    for ll in open(dirr+'/variant_reads'):
        if '#' not in ll:
            dpp=int(ll.split('\t')[9].split(':')[2])
            chrr=ll.split('\t')[0]
            pos=ll.split('\t')[1]
           

            read = ll.split('\t')[13].split('/')[0]
            if read not in count:
                count[read]=1
            if (chrr+'*'+pos) not in pos_map_reads:
                pos_map_reads[(chrr+'*'+pos)]=[read]
            else:
                pos_map_reads[(chrr+'*'+pos)].append(read)

            if read not in reads_map_pos:
                reads_map_pos[read]=[chrr+'*'+pos]
            else:
                reads_map_pos[read].append(chrr+'*'+pos)


    for p in pos_map_reads:
         var_dp[p]=len(pos_map_reads[p])

    
    fas=open(dirr+'/test4.fa','w')
    for l in open(dirr+'/tmp.fa'):
        l=l.rstrip()
        if l[:4]=='@SRR' and l[1:] in count:
            flag=1
            print('>'+l[1:],sep='',file=fas)
        elif flag==1:
            print(l,file=fas)
            flag=0
        else:
            flag=0
    
    os.unlink(dirr+'/tmp2_1.bam')
    os.unlink(dirr+'/variant_reads')
    os.unlink(dirr+'/tmp.fa')

    
    os.system('blat '+genome_file+' '+dirr+'/test4.fa '+dirr+'/blat.out -out=blast8')
    
    os.unlink(dirr+'/test4.fa')
    

    for li in open(dirr+'/blat.out'):
        li=li.rstrip()
        '''
        chr_pos_label=li.split('\t')[0].split('@')[0]
        chr_label=li.split('\t')[0].split('@')[0].split('_')[0]
        pos_label=int(li.split('\t')[0].split('@')[0].split('_')[1])
        '''
        reads=li.split('\t')[0]
        cov=float(li.split('\t')[-1])
        chrr=li.split('\t')[1]
        start=int(li.split('\t')[8])
        end=int(li.split('\t')[9])
       
        if count[reads]==1:
            count[reads]+=1
            for label in reads_map_pos[reads]:
                
                if chrr==label.split('*')[0] and int(label.split('*')[1])>=start and int(label.split('*')[1])<=end:
                    pass
                    
                else:
                    var_dp[label]-=1
        elif count[reads]==2:
            count[reads]+=1
            if cov >95.0 :
                for label in reads_map_pos[reads]:
                    var_dp[label]-=1 
        else:
            continue 



    os.unlink(dirr+'/blat.out')


    os.system('samtools view -b -f 0x80 '+dirr+'/tmp2.bam'+' -o '+dirr+'/tmp2_2.bam')



    os.system('bedtools intersect -wb -a '+dirr+'/filter2_variants.vcf'+' -b '+dirr+'/tmp2_2.bam'+' > '+dirr+'/variant_reads')

    os.system('bedtools bamtofastq -i '+dirr+'/tmp2_2.bam -fq '+dirr+'/tmp.fa')


    reads_map_pos={}
    pos_map_reads={}
    count={}
    for ll in open(dirr+'/variant_reads'):
        if '#' not in ll:
            dpp=int(ll.split('\t')[9].split(':')[2])
            chrr=ll.split('\t')[0]
            pos=ll.split('\t')[1]
           

            read = ll.split('\t')[13].split('/')[0]
            if read not in count:
                count[read]=1
            if (chrr+'*'+pos) not in pos_map_reads:
                pos_map_reads[(chrr+'*'+pos)]=[read]
            else:
                pos_map_reads[(chrr+'*'+pos)].append(read)

            if read not in reads_map_pos:
                reads_map_pos[read]=[chrr+'*'+pos]
            else:
                reads_map_pos[read].append(chrr+'*'+pos)

    for p in pos_map_reads:
         if p in var_dp:
             var_dp[p]+=len(pos_map_reads[p])
         else:
             var_dp[p]=len(pos_map_reads[p])




    fas=open(dirr+'/test4.fa','w')
    for l in open(dirr+'/tmp.fa'):
        l=l.rstrip()
        if l[:4]=='@SRR' and l[1:] in count:
            flag=1
            print('>'+l[1:],sep='',file=fas)
        elif flag==1:
            print(l,file=fas)
            flag=0
        else:
            flag=0


    os.unlink(dirr+'/tmp2_2.bam')
    os.unlink(dirr+'/variant_reads')
    os.unlink(dirr+'/tmp.fa')

    os.system('blat '+genome_file+' '+dirr+'/test4.fa '+dirr+'/blat.out -out=blast8')

    os.unlink(dirr+'/test4.fa')

    for li in open(dirr+'/blat.out'):
        li=li.rstrip()
        '''
        chr_pos_label=li.split('\t')[0].split('@')[0]
        chr_label=li.split('\t')[0].split('@')[0].split('_')[0]
        pos_label=int(li.split('\t')[0].split('@')[0].split('_')[1])
        '''
        reads=li.split('\t')[0]
        cov=float(li.split('\t')[-1])
        chrr=li.split('\t')[1]
        start=int(li.split('\t')[8])
        end=int(li.split('\t')[9])
       
        if count[reads]==1:
            count[reads]+=1
            for label in reads_map_pos[reads]:
                
                if chrr==label.split('*')[0] and int(label.split('*')[1])>=start and int(label.split('*')[1])<=end:
                    pass
                    
                else:
                    var_dp[label]-=1
        elif count[reads]==2:
            count[reads]+=1
            if cov >95.0 :
                for label in reads_map_pos[reads]:
                    var_dp[label]-=1 
        else:
            continue 

    os.unlink(dirr+'/blat.out')
    os.unlink(dirr+'/tmp2.bam')


    final=open(dirr+'/rna_editing.vcf','w')   
    for l in open(dirr+'/filter2_variants.vcf'):
        l=l.rstrip()
        if '#' in l:
            print(l,file=final)
        else:
            chrr=l.split('\t')[0]
            pos=l.split('\t')[1]
            if var_dp[chrr+'*'+pos]>0:
                print(l,file=final)


    os.unlink(dirr+'/filter1_variants.vcf')
    os.unlink(dirr+'/filter1_variants.vcf.idx')
    os.unlink(dirr+'/splited.log')
    os.unlink(dirr+'/filter2_variants.vcf')
    os.unlink(dirr+'/tmp.bam')
    os.unlink(dirr+'/tmp.bai')
    os.unlink(dirr+'/test4.bam.bai')
    os.unlink(dirr+'/test.txt')

    os.unlink(tmp_out4)


               
            


     
   
    
    
    



    
    


