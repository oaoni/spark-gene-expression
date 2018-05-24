import numpy as np
import pandas as pd
import json
import requests
import os
import re
from io import StringIO
from io import BytesIO
import tarfile

#Generates a folder to store the data portal gene expression data if none exits
newpath = os.path.join(os.getcwd(),"data")
if not os.path.exists(newpath):
    os.makedirs(newpath)

class gdc_mirna:
    """
    Creates data objects that can query the gdc data portal for miRNA expression data and
    read the data as a pandas dataframe (gene x sample_id).
    """

    def __init__(self, name):
        #Initialize the type of cancer for the database query and the size of query
        self.name = name
        #initialize gene epression data matrix
        self.data = pd.DataFrame()
        #Initialize an empty http reponse
        self.response = ''
        #Folder to store saved data
        self.main_dir = os.path.join(os.getcwd(),"data","miRNA")
        self.query_dir = os.path.join(self.main_dir,self.name)
        #Initialize location of data csv file
        self.file = ''
        #Initialize index
        self.index = ''

    def data_query(self):
        """
        Performs a query of the NCI genomic portal given a type of cancer initialized with the class.
        Ex. Type: Hepatocellular Carcinoma - LIHC
        Input: self.name followed by no. of samples desired. Ex. LIHC10 returns gene expression for
        first 10 samples. If specific number not present, will return all samples in database.
        Output: Binary data file of compressed tar.gz file in memory
        """
        files_endpt = "https://api.gdc.cancer.gov/files"

        #Parse name for type of cancer and desired number of samples integer
        cancer = re.search(r'\D+', self.name).group(0)
        size = re.search(r'\d+', self.name)
        if size:
            size = size.group(0)
        else:
            size = 2000

        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content":
                    {
                        "field": "cases.project.project_id",
                        "value": ["TCGA-"+cancer]
                    }
                },
                {
                    "op": "in",
                    "content":
                    {
                        "field": "files.experimental_strategy",
                        "value": ["miRNA-Seq"]
                    }
                },
                {
                    "op": "in",
                    "content":
                    {
                        "field": "files.data_type",
                        "value": ["miRNA Expression Quantification"]
                    }
                }
            ]
        }

        # Here a GET is used, so the filter parameters should be passed as a JSON string.
        params = {
            "filters": json.dumps(filters),
            "fields": "file_id",
            "format": "JSON",
            "size": size  #Set to the first 10 files for development
        }

        response = requests.get(files_endpt, params = params)
        file_uuid_list = []

        # This step populates the download list with the file_ids from the previous query
        for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
            file_uuid_list.append(file_entry["file_id"])

        data_endpt = "https://api.gdc.cancer.gov/data"

        params = {"ids": file_uuid_list}
        #Acquire memory location of compressed data from the data portal
        self.response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type":"application/json"})

    def data_read(self):
        """
        Extracts, decodes, and generates a pandas dataframe from the queried data
        Input: self.response.content
        Output: self.data stores a pandas dataframe (mirna x sample id)
        """
        #Run query if server response is empty
        if not self.response:
            self.data_query()

        #Uncompression of the data response
        with BytesIO(self.response.content) as targz:
            #Open a tarfile with response content
            with tarfile.open(fileobj=targz) as tar:
                for member in tar.getmembers()[1:]: #Iterate through the members of the tarfile
                    #extract,decode and concatenate member to the dataframe
                    with StringIO(tar.extractfile(member).read().decode('utf-8')) as data:
                        self.data = pd.concat([self.data,pd.read_table(data,sep="\t",usecols=['read_count']).
                                        rename(columns={'read_count':member.name.split('/')[0]})],axis=1)

                #index = pd.read_table(StringIO(tar.extractfile(tar.getmembers()[1]).read().decode('utf-8')),
               # sep="\t",usecols=['miRNA_ID'])
                self.data.index = pd.read_table(StringIO(tar.extractfile(tar.getmembers()[1]).
                                                         read().decode('utf-8')),sep="\t",
                                                usecols=['miRNA_ID']).miRNA_ID.tolist()
                self.data.index.name = 'miRNA_ID'

    def data_save(self, safe=True, format="csv"):
        """
        Saves loaded data as a csv, txt or in parquet format
        Inputs: safe = True/False, format = ['csv','txt']
        Output: Saved file in respective folder in the query_dir
        """
        #Create a path to save the data if it doesnt exist already
        if not os.path.exists(self.query_dir):
            os.makedirs(self.query_dir)

        if self.data.empty:
            print('Data has not been queried or read, run method self.data_read')
            return
        elif format == "csv":
            self.file = os.path.join(self.query_dir,self.name+".csv")
            self.data.to_csv(self.file)
            print("csv file successfully saved...")
        elif format == "txt":
            self.file = os.path.join(self.query_dir,self.name+".txt")
            self.data.to_csv(self.file,sep='\t')
            print("txt file successfully saved...")
        # elif format == "parquet":
        #     self.file = os.path.join(self.query_dir,self.name+".pq")
        #     self.data.to_parquet(self.file)
        #     print("parquet file successfully saved")
        # elif format == "sql":
        #     self.file = os.path.join(self.query_dir,self.name)
        #     self.data.to_sql(self.file)
        #     print("sql file successfully saved")

if __name__ == '__main__':


    KIRC = gdc_mirna('KIRC')
    KIRC.data_read()
    print(KIRC.data.shape)
    KIRC.data_save(format="csv")
