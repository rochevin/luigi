import json
import subprocess
import os

###Usefull fonctions for each pipeline

##Launch a bash command in python, return the value of the command, for instance ls return a string with all files in current directory
def run_cmd(cmd):
	p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
	ret_code = p.wait()
	output = p.communicate()[0]
	return output

##Open a json file and return a dictionnary with all values
def json_read(json_file):
	json_data = open(json_file)
	return json.load(json_data)

##return true if the given file name exist
def exist(file_name):
	return os.path.isfile(file_name)


###All NGS python class

##class with all command for samtools
##dependencies : 
class SamtoolsCommand(object):
	"""All command from samtools inside methods :
	-samtools_index_ref
	-samtools_index_bam
	-samtools_sam_to_bam
	-samtools_sort_bam
	-samtools_remove_PCR_duplicate
	-is_indexing : check if the given file is indexed given an extension"""
	def samtools_index_ref(self,fasta):
		"""Use samtools with samtools faidx"""
		return run_cmd(["samtools","faidx",fasta]) if not self.is_indexing(fasta) else None

	def samtools_sam_to_bam(self,file_index,file_output,file_input,options = ["-b","-S","-q","25"]):
		"""Use samtools with samtools view"""
		return run_cmd(["samtools","view"]+options+["-t",file_index,"-o",file_output,file_input]) if not exist(file_output) else None

	def samtools_sort_bam(self,bam_file,sorted_bam_file):
		"""Use samtools with samtools sort"""
		#check if the output file is without .bam extension (because samtools add the .bam automaticly for the output)
		file_name, extension = os.path.splitext(sorted_bam_file)
		good_output = sorted_bam_file if not extension else file_name
		return run_cmd(["samtools","sort",bam_file,good_output]) if not exist(file_name+extension) else None
	
	def samtools_remove_PCR_duplicate(self,bam_file,nodups_bam_file,opt="-s"):
		"""Use samtools with samtools rmdups"""
		return run_cmd(["samtools","rmdup",opt,bam_file,nodups_bam_file]) if not exist(nodups_bam_file) else None

	def samtools_index_bam(self,bam_file,index_bam_file):
		"""Use samtools with samtools index"""
		return run_cmd(["samtools","index",bam_file,index_bam_file]) if not exist(index_bam_file) else None

	def samtools_merge_bam(self,merged_bam_file,bam_list):
		"""use samtools with samtools merge"""
		return run_cmd(["samtools","merge",merged_bam_file]+list(bam_list))


	def is_indexing(self,initial_file,extensions=[".fai"]):
		"""Check if the ref genome is indexed, defaut : [".fai"]"""
		result_bool = [os.path.isfile(initial_file+extension) for extension in extensions ]
		return True in result_bool

	def bam_for_sample(self,sample_list):
		return [sample+".bam" for sample in sample_list]

class BwaCommand(object):
	""""""
	def bwa_index(self,genome):
		"""Use BWA with bwa index to index fasta file"""
		return run_cmd(["bwa","index",genome]) if not self.is_indexing(genome) else None


	def bwa_aln(self,fastq_list,output,genome,threads = "2"):
		"""Use BWA with bwa aln for fastq file(s)"""
		fastq_list = self.get_fastq_for_sample([fastq_list]) if type(fastq_list) == "str" else self.get_fastq_for_sample(list(fastq_list))
		command = ["bwa",
						"aln",
						"-t",
						threads,
						"-f",
						output,
						genome]+fastq_list
		print(command)
		return run_cmd(command) if not exist(output) else None

	def bwa_samse(self,fastq_list,genome,file_input,file_output):
		"""Use BWA with bwa samse for fastq file(s)"""
		fastq_list = self.get_fastq_for_sample([fastq_list]) if type(fastq_list) == "str" else self.get_fastq_for_sample(list(fastq_list))
		return run_cmd(["bwa",
						 "samse",
						 "-f",
						 file_output,
						 genome,
						 file_input
						 ]+fastq_list) if not exist(file_output) else None

	def bwa_sampe(self,sai_list,fastq_list,genome,file_output):
		"""Use BWA with bwa sampe for 2 fastq files"""
		sai_list = self.get_sai_for_sample([sai_list]) if type(sai_list) == "str" else self.get_sai_for_sample(list(sai_list))
		fastq_list = self.get_fastq_for_sample([fastq_list]) if type(fastq_list) == "str" else self.get_fastq_for_sample(list(fastq_list))
		return run_cmd(["bwa",
						 "sampe",
						 "-f",
						 file_output,
						 genome
						 ]+sai_list+fastq_list) if not exist(file_output) else None


	def is_indexing(self,initial_file,extensions=[".amb",".ann",".bwt",".pac",".sa"]):
		"""Check if the ref genome is indexed, defaut : [".amb",".ann",".bwt",".pac",".sa"]"""
		result_bool = [os.path.isfile(initial_file+extension) for extension in extensions ]
		return True in result_bool

	def get_fastq_for_sample(self,sample_list):
		return [sample+".fastq" for sample in sample_list]

	def get_sai_for_sample(self,sample_list):
		return [sample+".sai" for sample in sample_list]


class Cutadapt(object):
	"""Use cutadapt with -g/-a option for 5' or 3' adapter"""
	def cutadapt_command(self,adapter,file_output,file_input,opt = "-g"):
		return run_cmd(["cutadapt",
						opt,
						 adapter,
						 "-o",
						 file_output,
						 file_input
						 ])
