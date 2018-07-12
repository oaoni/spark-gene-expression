import numpy as np
import pandas as pd
import requests
import os
import re
import json
from io import StringIO
from io import BytesIO
import tarfile
import sqlite3

class gdc_cnv:
    """
    Creates data objects that can query the gdc data portal for copy number variation data
    """

    def __init__(self, name):
        #Initialize the type of cancer for the database query and the size of query
        self.name = name
        #Initialize an empty http reponse
        self.response = ''
        #Folder to store saved data
        self.main_dir = os.path.join(os.getcwd(),"data","CNV")
        self.query_dir = os.path.join(self.main_dir,self.name)
        #Initialize location of data csv file
        self.file = ''
        #Initialize index
        self.index = ''

    def data_query(self):
        """
        Performs a query of the NCI genomic portal given a type of cancer initialized with the class.
        Ex. Type: Hepatocellular Carcinoma - LIHC
        Input: self.name followed by no. of samples desired. Ex. LIHC10 returns cnv data for
        first 10 samples. If specific number not present, will return all samples in database.
        Output: Binary data file of compressed tar.gz file in memory as self.response
        """

        files_endpt = "https://api.gdc.cancer.gov/files"

        #Parse name for type of cancer and desired number of samples integer
        cancer = re.search(r'\D+', self.name).group(0)
        size = re.search(r'\d+', self.name)
        if size:
            size = size.group(0)
        else:
            size = 2000

        #Specify the type of copy number segmentation for the query
        CNS = "Copy Number Segment" #Or ["Masked Copy Number Segment"]

        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "cases.project.project_id",
                        "value": ["TCGA-"+cancer]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "files.data_type",
                        "value": [CNS]
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
        Extracts, decodes, and generates a
        Input: self.response.content
        Output:
        """
        #Run query if server response is empty
        if not self.response:
            self.data_query()

        #Connect to the sqlite database
        connection = sqlite3.connect(os.path.join(self.main_dir,self.name+'_cnv.sqlite'))

        #Access body of bytes response and open compressed targz file as a tarfile
        with tarfile.open(fileobj=BytesIO(self.response.content)) as targz:
            #Iterate through the members of the tarfile
            for member in targz.getmembers()[1:]:
                #Extract each member
                with BytesIO(targz.extractfile(member).read()) as data:
                    #Read member into a pandas dataframe and store in a sql database
                    pd.read_table(data, sep='\t').to_sql(
                        name = member.name.split('/')[0],
                        con = connection,
                        schema = 'GDC_Aliquot TEXT, Chromosome TEXT, Start INTEGER, End INTEGER, Num_Probes INTEGER, Segment_Mean REAL',
                        index=False,
                        if_exists='append'
                    )

if __name__ == '__main__':

    #Make the query object, while initializing the cancer type and the number of examples
    cnvKIRP = gdc_cnv('KIRP5')

    #Run query method to make https request to the data portal
    cnvKIRP.data_query()

    #Read the contents of the query
    cnvKIRP.data_read()
