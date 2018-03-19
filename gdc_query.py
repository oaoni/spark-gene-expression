import requests
import json
import re
import gzip
import pandas as pd
import tarfile
import os
import sys
import shutil
import numpy as np
import tarfile
import io
from io import StringIO

#Generates a folder to store the data portal gene expression data if none exits
newpath = os.path.join(os.getcwd(),"data")
if not os.path.exists(newpath):
    os.makedirs(newpath)

class gdc_data:
    '''
    Creates data objects that can query the gdc data portal for gene expression data,
    write compressed data from portal to disk, and uncompress and store gene expression
    data in a pandas dataframe (gene x sample_id)
    '''

    def __init__(self, name):
        #Initialize the type of cancer for the database query and the size of query
        self.name = name
        #initialize gene epression data matrix
        self.data = pd.DataFrame()
        #Initialize empty manifest data matrix
        self.manifest = pd.DataFrame()
        #Initialize the location for the data directory
        self.main_dir = os.path.join(os.getcwd(),"data")
        self.query_dir = os.path.join(self.main_dir,self.name)
        #Initialize location of data csv file
        self.file = os.path.join(self.query_dir,self.name+".csv")
        #Initialize empty file name
        self.file_name = ''
        #Initialize an empty http reponse
        self.response = ''
        #Initialize variable for size of query
        self.size = ''

    def data_query(self):
        '''
        Performs a query of the NCI genomic portal given a type of cancer.
        Ex. Type: Hepatocellular Carcinoma - LIHC
        Name followed by no. of samples desired. Ex. LIHC10 returns gene expression for first 10 samples
        Returns the name of compressed tar.gz file, and a binary data file in memory
        '''

        files_endpt = "https://api.gdc.cancer.gov/files"

        cancer = ''.join([x for x in self.name if not x.isdigit()])
        size = ''.join([x for x in self.name if x.isdigit()])
        #If the size of the query was not specified, acquire data for all samples
        if not size.isdigit():
            size = 2000

        #Filters for the query, recieving all RNA-Seq, HTSeq-Count files for a specific cancer
        filters = {
            "op": "and",
            "content":[
                {
                "op": "in",
                "content":{
                    "field": "cases.project.project_id",
                    "value": ["TCGA-"+cancer]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq"]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - Counts"]
                    }
                }
            ]
        }

        # Here a GET is used, so the filter parameters should be passed as a JSON string.
        params = {
            "filters": json.dumps(filters),
            "fields": "file_id",
            "format": "JSON",
            "size": size  #Set to the first 10 files for developing
            }

        response = requests.get(files_endpt, params = params)
        file_uuid_list = []

        # This step populates the download list with the file_ids from the previous query
        for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
            file_uuid_list.append(file_entry["file_id"])

        data_endpt = "https://api.gdc.cancer.gov/data"

        params = {"ids": file_uuid_list}
        #Acquire memory location of compressed data from the data portal
        response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

        response_head_cd = response.headers["Content-Disposition"]
        #Acquire the name of the file
        file_name = re.findall("filename=(.+)", response_head_cd)[0]

        self.file_name = file_name
        self.response = response

    def data_write(self):
        """
        Writes query file to query directory
        """

        #Performs data query if filename and response have not been populated yet
        if not self.file_name and not self.response:
            self.data_query()

        #Create a path for this query if it doesnt exist already
        if not os.path.exists(self.query_dir):
            os.makedirs(self.query_dir)

        #desired location of compressed data targz file
        targz = os.path.join(self.query_dir,self.file_name)
        #Opens a file named after the file_name, and writes the contents of the query to the file
        with open(targz, "wb") as output_file:
            output_file.write(self.response.content) #writes the response to the desired location

        self.data_write_targz()

    def data_write_targz(self):
        '''
        Uncompressess a targz file under the query directory, and writes to disk
        '''

        #Performs data query and writes targz if filename and response have not been populated yet
        #Check if the filename has been populated, if not then check if it exists in the query directory
        #, if neither, run query
        if not self.file_name:
            self.file_name = ''.join([x for x in os.listdir(self.query_dir) if x[-6:] == 'tar.gz'])
            if not self.file_name:
                self.data_query()
                self.data_write()

        #desired location of compressed data targz file
        targz = os.path.join(self.query_dir,self.file_name)

        #Create a path for the uncompressed targz files if it doesnt exist already
        uncomp_targz_dir = os.path.join(self.query_dir,"uncompressed_targz")
        if not os.path.exists(uncomp_targz_dir):
            os.makedirs(uncomp_targz_dir)
        #Unzips the tar.gz file into desired folder
        with tarfile.open(targz) as tar:
            tar.extractall(uncomp_targz_dir)
            tar.close()
        #Stores the manifest of the data
        self.manifest = pd.read_table(os.path.join(uncomp_targz_dir,"MANIFEST.txt"),sep="\t")

        #Create a path for this uncompressed gz files if it doesnt exist already
        uncomp_gz_dir = os.path.join(uncomp_targz_dir,"uncompressed_gz")
        if not os.path.exists(uncomp_gz_dir):
            os.makedirs(uncomp_gz_dir)

        #Unzips all gz gene expression files in the query directory
        for subdir, dirs, files in os.walk(self.query_dir):
            for file in files:
                if file[-4:] == "s.gz":
                    with gzip.open(os.path.join(subdir,file),'rb') as f:
                        file_content = f.read().decode("utf-8")
                        df = pd.read_csv(StringIO(file_content),sep="\t",header=None).set_index(0)
                        df.columns = [files[0]]
                        df.to_csv(os.path.join(uncomp_gz_dir,file[:-3]),header=False,sep=",",index=True)

        self.data_save()

    def data_save(self):
        """
        Saves dataframe to object from uncompressed tar and targz files
        """
        #initialize/clear gene epression data matrix
        self.data = pd.DataFrame()

        uncomp_targz_dir = os.path.join(self.query_dir,"uncompressed_targz")
        uncomp_gz_dir = os.path.join(uncomp_targz_dir,"uncompressed_gz")
        #Stores the manifest of the data
        self.manifest = pd.read_table(os.path.join(uncomp_targz_dir,"MANIFEST.txt"),sep="\t")

        for subdir, dirs, files in os.walk(self.query_dir):
            for file in files:
                if file[-4:] == "unts":
                    df = pd.read_csv(os.path.join(uncomp_gz_dir,file),sep=",",header=None).set_index(0)
                    df.columns = [file]
                    self.data = pd.concat([self.data,df],axis=1)

        self.size = self.data.shape #Store the dimensions of the data matrix

    def data_add(self,data):
        """
        Directly load data into object if available as a dataframe
        """
        self.data = data

    def save_file(self,safe=True):
        """
        Saves loaded data as a csv to disk
        safe mode (default) - True
        if downloading and saving large volumes of data change to (False, 0, None,"",[])
        """
        #File location
        file = os.path.join(self.query_dir,self.name+".csv")

        if safe:
            #Overwrite prevention
            if os.path.exists(file):
                response = input("Do you want to overwrite previously saved file? (y/n)")
                if response == "y":
                    self.data.to_csv(file)
                    print("csv file successfully saved")
            else:
                self.data.to_csv(file)
                print("csv file successfully saved")

            input("Press Enter to Continue...")
        else:
            #No overwrite protection, but no user prompts for loops
            self.data.to_csv(file)
            print("csv file successfully saved")

    def read_csv(self):
        """
        Reads data from csv to pandas dataframe
        """

        if os.path.exists(self.file):
            self.data = pd.read_csv(self.file)
        else:
            print("file does not exist")
			
#Debugging and testing		
if __name__ == "__main__":

	#Initialize a data object
	test = gdc_data("LIHC10")
	
	#Read data from a csv file in query directory
	test.read_csv()
	
	#Head of loaded data
	print(test.data.head())
	