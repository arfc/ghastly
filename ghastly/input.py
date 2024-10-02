import json

class Input:
    '''
    class for containing data from a ghastly input file (json)
    '''
    def __init__(self, input_file):
        '''
        initializes ghastly input object using a ghastly input file.
        '''

        f = open(input_file, 'r')
        params = json.load(f)
        
        self.simulation = params["simulation"]
        self.core_intake = params["core_intake"]
        self.core_main = params["core_main"]
        self.core_outtake = params["core_outtake"]

