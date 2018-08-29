import numpy as np
import pandas as pd
import requests
import os
import re
import json
from io import StringIO
from io import BytesIO
import tarfile
import gzip

class gdc_snv:
    """
    Creates data objects that can query the gdc data portal for Simple Nucleotide
    Variation data
    """

    def __init__(self, name):
        #Initialize the type of cancer for the database query and the size of query
        self.name = name
        #Initialize an empty http reponse
        self.response = ''
        #Folder to store saved data
        self.main_dir = os.path.join(os.getcwd(),"data","SNV")
        #Initialize variable to store query manifest
        self.manifest = pd.DataFrame()
        #Initialize location of data csv file
        self.file = ''
        #Initialize variable for dataframe
        self.data = ''

    def data_query(self):
        """
        Performs a query of the NCI genomic portal given a type of cancer initialized with the class.
        Ex. Type: Hepatocellular Carcinoma - LIHC
        Input: self.name
        Output: Binary data file of compressed tar.gz file in memory as self.response
        """

        files_endpt = "https://api.gdc.cancer.gov/files"

        #Generate filters for https request
        filters = {
            "op":"and",
            "content":[
                {
                    "op":"in",
                    "content":{
                        "field": "cases.project.project_id",
                        "value": ["TCGA-"+self.name]
                    }
                },
                {
                    "op":"in",
                    "content":{
                        "field": "files.analysis.workflow_type",
                        "value": ["MuTect2 Variant Aggregation and Masking"]
                    }
                },
                {
                    "op":"in",
                    "content":{
                        "field": "files.data_type",
                        "value": ["Masked Somatic Mutation"]
                    }
                }
            ]
        }

        # Here a GET is used, so the filter parameters should be passed as a JSON string.
        params = {
            "filters": json.dumps(filters),
            "fields": "file_id",
            "format": "JSON",
            "size": 1  #Acquires only the first file
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
        Output: self.data stores a pandas dataframe
        """

        #Run query if server response is empty
        if not self.response:
            self.data_query()


        with gzip.open(BytesIO(self.response.content), 'rb') as gz:
            with StringIO(gz.read().decode('utf-8')) as data:
                self.data = pd.read_csv(data, sep="\t", skiprows=range(0,5),
                usecols = list(range(0,88))+list(range(90,98))+list(range(99,120)),
                low_memory = False)

    def data_save(self, format="csv"):
        """
        Saves loaded data as a csv or txt file
        Inputs: format = ['csv','txt']
        Output: Saved file in respective folder in the query_dir
        """

        if self.data.empty:
            print('Data has not been queried or read, run method self.data_read')
            return
        elif format == "csv":
            self.file = os.path.join(self.main_dir,self.name+"_SNV.csv")
            self.data.to_csv(self.file)
            print("csv file successfully saved...")
        elif format == "txt":
            self.file = os.path.join(self.main_dir,self.name+"_SNV.txt")
            self.data.to_csv(self.file,sep='\t')
            print("txt file successfully saved...")


if __name__ == "__main__":

    #kirc = gdc_snv("KIRC")
    #kirc.data_read()
    #print(kirc.data.shape)
    #kirc.data_save()

    names = ['KIRP','LIHC','LUAD','LUSC','CHOL','PRAD','THCA','HNSC']

    for name in names:
        type = gdc_snv(name)
        type.data_read()
        type.data_save()
