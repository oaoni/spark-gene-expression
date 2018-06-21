import requests
import json
import re
import gzip
import pandas as pd
import tarfile
import os
from io import StringIO
from io import BytesIO
import time


#Generates a folder to store the data portal gene expression data if none exits
newpath = os.path.join(os.getcwd(),"data")
if not os.path.exists(newpath):
    os.makedirs(newpath)

class gdc_rnaseq:
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
        self.main_dir = os.path.join(os.getcwd(),"data","RNASeq")
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
        response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

        response_head_cd = response.headers["Content-Disposition"]
        #Acquire the name of the file
        file_name = re.findall("filename=(.+)", response_head_cd)[0]

        self.file_name = file_name
        self.response = response

    def data_read(self):
        """
        Extracts, decodes, and generates a pandas dataframe from the queried data
        Input: self.response.content
        Output: self.data stores a pandas dataframe (mirna x sample id)
        """
        #Run query if server response is empty
        if not self.response:
            self.data_query()

        #Store the contents of the query in a memory location
        with BytesIO(self.response.content) as targz:
            #Generate a tarfile from the compressed bytes object
            with tarfile.open(fileobj=targz) as tar:
                #Iterate through the members of the targz file
                for member in tar.getmembers()[1:]:
                    #Extract member of the targz file, read the data, and store the gz file in bytes memory location,
                    #and open gz file
                    with gzip.open(BytesIO(tar.extractfile(member).read())) as gz:
                        #Read gz file, decode to utf-8 and set in string memory location
                        with StringIO(gz.read().decode("utf-8")) as data:
                            #Read the txt data file and concatenate the the dataframe
                            self.data = pd.concat([self.data,pd.read_table(data,sep="\t", header=None,usecols=[1])
                            .rename(columns={1:member.name.split('/')[0]})],axis=1)

                #Set index of rnaseq names on the dataframe
                self.data.index = pd.read_table(StringIO(gzip.open(BytesIO((tar.extractfile(tar.getmembers()[1])
                .read()))).read().decode("utf-8")),sep="\t",header=None,usecols=[0])[0].tolist()
                #Set index name
                self.data.index.name = 'RNASeq_ID'

    def data_save(self, safe=True, format="csv"):
        """
        Saves loaded data as a csv, txt or in parquet format
        Inputs: safe = True/False, format = ['csv','txt']
        Output: Saved file in respective folder in the query_dir
        """
        #Create a path to save the data if it doesnt exist already
        if not os.path.exists(self.main_dir):
            os.makedirs(self.main_dir)

        if self.data.empty:
            print('Data has not been queried or read, run method self.data_read')
            return
        elif format == "csv":
            self.file = os.path.join(self.main_dir,self.name+"_RNASeq.csv")
            self.data.to_csv(self.file)
            print("csv file successfully saved...")
        elif format == "txt":
            self.file = os.path.join(self.main_dir,self.name+"_RNASeq.txt")
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

        self.data_store()

    def data_store(self):
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
            self.data = pd.read_csv(self.file,index_col=0)
        else:
            print("file does not exist")

#Debugging and testing
if __name__ == "__main__":

    t0 = time.time()
    #Initialize a data object
    KIRP = gdc_rnaseq("KIRC")
    KIRP.data_read()
    KIRP.data_save()
    print(KIRP.data.shape)
    t1 = time.time()

    print(t1-t0)

    #print(test.data.shape)
    #test.data_read()
