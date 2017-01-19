import subprocess
import shlex

#PATHS_path = '/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/variantannotation/variantannotation/PATHS.txt'
def run_annovar(annovar_path, input_vcf_path, output_csv_path):

#args_str = "sudo perl /database/annovar/table_annovar.pl /data/Nof1/normal_blood_WGS.vqsr.vcf /database/annovar/humandb/ -buildver hg19 -out /data/ccbb_internal/interns/Carlo/annovar_out/SUBPROCESS -remove -protocol knownGene,tfbsConsSites,cytoBand,targetScanS,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all,cg46,cg69,popfreq_all,clinvar_20140929,cosmic70,nci60 -operation g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -csvout"

    args_str = "sudo perl " + annovar_path + "table_annovar.pl " + input_vcf_path + " " + annovar_path + "humandb/ -buildver hg19 -out " + output_csv_path + " -remove -protocol knownGene,tfbsConsSites,cytoBand,targetScanS,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all,cg46,cg69,popfreq_all,clinvar_20140929,cosmic70,nci60 -operation g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -csvout"
    args = shlex.split(args_str)
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    return p.communicate()




