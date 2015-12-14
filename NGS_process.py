import json
import subprocess
import os

###Usefull fonctions for each pipeline
##Parse .conf file for NGS pipeline

def read_conf(file_name,main_opts=["GENOME","FASTQ","DIRECTORY"],sub_opt=["NAME","DIR"]):

	config = {"DIRECTORY":[],"NAMES":[]}
	content = []
	global_option = ""
	aln_option = ""
	name = ""
	directory = ""
	with open(file_name,"r") as conf_file:
		for line in conf_file:
			line = line.replace("\t","")
			if line.startswith("#") or line.startswith("\n"):
				continue
			else:
				line = line.replace("\n","")
				if line.startswith("["):
					if content and global_option in main_opts:
						if global_option == "FASTQ" and name:
							config[name] = content
							content = [] if option == "NAME" else content
						elif global_option == "GENOME":
							config[global_option] = content
							content = []
						elif global_option == "DIRECTORY":
							config[global_option] = content
							content = []
						
					result = line.replace("[","").replace("]","").split(":")
					option = result[0].upper()
					if option in main_opts:
						global_option = option
					elif option in sub_opt:
						if option == "NAME":
							name = result[1]
							config["NAMES"].append(name)
						elif option == "DIR":
							directory = result[1]
				else:
					if global_option:
						if global_option == "FASTQ" and name:
							result = line.replace("\t","").replace(" ","").replace("\"","").split(",")
							content.append([(directory,fichier) for fichier in result if fichier])
						elif global_option in ["GENOME","DIRECTORY"] :
							result = line.replace("\t","").replace(" ","").replace("\"","").split(":")
							try:
								content.append(result[1])
							except IndexError:
								content.append(result[0])
		if content:
			if global_option == "FASTQ" and name:
				config[name] = content
			elif global_option == "GENOME":
				config[global_option] = content
			elif global_option == "DIRECTORY":
				config[global_option] = content

				
		for key,first_dim in config.items():
			if key in config["NAMES"]:
				res = []
				for sec_dim in first_dim:
					sub_res = []
					for third_dim in sec_dim:
						directory = int(third_dim[0])-1
						file_name = third_dim[1]
						try:
							path = config["DIRECTORY"][directory]
						except IndexError:
							path = ""
						sub_res.append(os.path.join(path,file_name))

					res.append(sub_res)
				config[key] = res
		conf_file.close()
	return config




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


##add_directory to list
def add_directory(directory,liste):
	return [os.path.join(directory,elmt) for elmt in liste]

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
		command = ["samtools","faidx",fasta]
		print(" ".join(command))
		return run_cmd(command) if not self.is_indexing(fasta) else None

	def samtools_sam_to_bam(self,file_index,file_output,file_input,options = ["-b","-S","-q","25"]):
		"""Use samtools with samtools view"""
		command = ["samtools","view"]+options+["-t",file_index,"-o",file_output,file_input]
		print(" ".join(command))
		return run_cmd(command) if not exist(file_output) else None

	def samtools_sort_bam(self,bam_file,sorted_bam_file):
		"""Use samtools with samtools sort"""
		#check if the output file is without .bam extension (because samtools add the .bam automaticly for the output)
		file_name, extension = os.path.splitext(sorted_bam_file)
		good_output = sorted_bam_file if not extension else file_name
		command = ["samtools","sort",bam_file,good_output]
		print(" ".join(command))
		return run_cmd(command) if not exist(file_name+extension) else None
	
	def samtools_remove_PCR_duplicate(self,bam_file,nodups_bam_file,opt="-s"):
		"""Use samtools with samtools rmdups"""
		command = ["samtools","rmdup",opt,bam_file,nodups_bam_file]
		print(" ".join(command))
		return run_cmd() if not exist(nodups_bam_file) else None

	def samtools_index_bam(self,bam_file,index_bam_file):
		"""Use samtools with samtools index"""
		command = ["samtools","index",bam_file,index_bam_file]
		print(" ".join(command))
		return run_cmd(command) if not exist(index_bam_file) else None

	def samtools_merge_bam(self,merged_bam_file,bam_list):
		"""use samtools with samtools merge"""
		command = ["samtools","merge",merged_bam_file]+list(bam_list)
		print(" ".join(command))
		return run_cmd(command)


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
		command = ["bwa","index",genome]
		print(" ".join(command))
		return run_cmd() if not self.is_indexing(genome) else None


	def bwa_aln(self,fastq_list,output,genome,threads = "2"):
		"""Use BWA with bwa aln for fastq file(s)"""
		fastq_list = self.get_fastq_for_sample([fastq_list]) if type(fastq_list) == str else self.get_fastq_for_sample(list(fastq_list))
		command = ["bwa",
						"aln",
						"-t",
						threads,
						"-f",
						output,
						genome]+fastq_list
		print(" ".join(command))
		return run_cmd(command) if not exist(output) else None

	def bwa_samse(self,fastq_list,genome,file_input,file_output):
		"""Use BWA with bwa samse for fastq file(s)"""
		fastq_list = self.get_fastq_for_sample([fastq_list]) if type(fastq_list) == str else self.get_fastq_for_sample(list(fastq_list))
		command = ["bwa",
						 "samse",
						 "-f",
						 file_output,
						 genome,
						 file_input
						 ]+fastq_list
		print(" ".join(command))
		return run_cmd(command) if not exist(file_output) else None

	def bwa_sampe(self,sai_list,fastq_list,genome,file_output):
		"""Use BWA with bwa sampe for 2 fastq files"""
		sai_list = self.get_sai_for_sample([sai_list]) if type(sai_list) == str else self.get_sai_for_sample(list(sai_list))
		fastq_list = self.get_fastq_for_sample([fastq_list]) if type(fastq_list) == str else self.get_fastq_for_sample(list(fastq_list))
		command = ["bwa",
						 "sampe",
						 "-f",
						 file_output,
						 genome
						 ]+sai_list+fastq_list
		print(" ".join(command))
		return run_cmd(command) if not exist(file_output) else None


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
