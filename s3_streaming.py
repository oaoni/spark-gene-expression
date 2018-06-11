import numpy as np
import pandas as pd
from io import StringIO
from io import BytesIO
import boto3

def numpy_to_s3(client,data,bucket,key):
    #Create a csv stream
    csv_buffer = StringIO()
    #Save numpy array as a csv to stream
    np.savetxt(csv_buffer, data, delimiter=',')
    #Convert the str to bytes
    body = bytes(csv_buffer.getvalue())
    #Save body to s3 in cooresponding bucket and key
    client.put_object(Body=body, Bucket=bucket, Key=key)
    
def s3_to_numpy(client,bucket,key):
    #Acquire object from s3
    obj = client.get_object(Bucket=bucket, Key=key)
    #Store read data from object in a bytes stream
    bytes_stream = BytesIO(obj['Body'].read())
    #Load data into a numpy array
    np_array = np.loadtxt(bytes_stream, delimiter=',')
    
    return np_array
	
def s3_to_pandas(client ,key, bucket):
    obj = client.get_object(Key=key, Bucket=bucket)
    bytes_stream = BytesIO(obj['Body'].read())
    data = pd.read_csv(bytes_stream, index_col=0)
    return data