'''
Created on Jul 8, 2014

@author: niels
'''



from thrift.transport import TSocket
from thrift.transport import TTransport
from thrift.protocol import TBinaryProtocol
from thrift.server import TServer

from openda_bmi.openda.bmi.thrift.BMIService import Iface
from openda_bmi.openda.bmi.thrift.BMIService import Processor
from openda_bmi.openda.bmi.thrift.ttypes import ModelException, BmiGridType

import sys
import signal
import numpy as np
import logging
import importlib

logger = logging.getLogger(__name__)

class ModelHandler(Iface):

    def __init__(self, model):
        self.model = model

    def initialize(self, config_file):
        if config_file is None or config_file == "":
            config_file = None
        
        try:
            model.initialize(config_file)
        except Exception as e:
            raise ModelException(str(e)) 
    
    def update(self):
        try:
            model.update()
        except Exception as e:
            raise ModelException(str(e)) 
        
    def update_until(self, time):
        try:
            model.update_until(time)
        except Exception as e:
            raise ModelException(str(e)) 


    def update_frac(self, frac):
        try:
            model.update_frac(frac)
        except Exception as e:
            raise ModelException(str(e))            
     
    def finalize_model(self):
        try:
            model.finalize()
        except Exception as e:
            raise ModelException(str(e)) 
     
    def get_component_name(self):
        try:
            return model.get_component_name()
        except Exception as e:
            raise ModelException(str(e)) 
    
    def get_input_var_names(self):
        try:
            return model.get_input_var_names()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_output_var_names(self):
        try:
            return model.get_output_var_names()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_var_type(self, long_var_name):
        """
        Get type of a variable, as a numpy dtype value string
        """
        try:
            return model.get_var_type(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_var_units(self, long_var_name):
        try:
            return model.get_var_units(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_var_rank(self, long_var_name):
        try:
            return model.get_var_rank(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_var_size(self, long_var_name):
        try:
            return model.get_var_size(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_var_nbytes(self, long_var_name):
        try:
            return model.get_var_nbytes(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_start_time(self):
        try:
            return model.get_start_time()
        except Exception as e:
            raise ModelException(str(e))

    def get_current_time(self):
        try:
            return model.get_current_time()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_end_time(self):
        try:
            return model.get_end_time()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_time_step(self):
        try:
            return model.get_time_step()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_time_units(self):
        try:
            return model.get_time_units()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_value(self, long_var_name):
        try:
            return model.get_value(long_var_name).flatten().tostring()
        except Exception as e:
            raise ModelException(str(e))
    
    def get_value_at_indices(self, long_var_name, inds):
        try:
            return model.get_value_at_indices(long_var_name, inds).tostring()
        except Exception as e:
            raise ModelException(str(e))
    
    def set_value(self, long_var_name, src):
        try:
            logger.info("received %s bytes", len(src))
            
            vartype = model.get_var_type(long_var_name)
            varshape = model.get_grid_shape(long_var_name)
            
            logger.info("var type, shape %s, %s", str(vartype), str(varshape))
            
            flatarray = np.fromstring(src, dtype=np.dtype(vartype))
            
            logger.info("flat array now shaped %s", str(flatarray.shape))
    
            value = np.reshape(flatarray, varshape)
            
            model.set_value(long_var_name, value)
        except Exception as e:
            raise ModelException(str(e))
    
    def set_value_at_indices(self, long_var_name, inds, src):
        try:
            vartype = model.get_var_type(long_var_name)
        
            model.set_value_at_indices(long_var_name, inds, np.fromstring(src, dtype=vartype))
        except Exception as e:
            raise ModelException(str(e))

    def get_grid_type(self, long_var_name):
        try:
            return model.get_grid_type(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_grid_shape(self, long_var_name):
        try:
            return model.get_grid_shape(long_var_name)
        except Exception as e:
            raise ModelException(str(e))

    def get_grid_spacing(self, long_var_name):
        try:
            return model.get_grid_spacing(long_var_name)
        except Exception as e:
            raise ModelException(str(e))

    def get_grid_origin(self, long_var_name):
        try:
            return model.get_grid_origin(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_grid_x(self, long_var_name):
        try:
            return model.get_grid_x(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_grid_y(self, long_var_name):
        try:
            return model.get_grid_y(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_grid_z(self, long_var_name):
        try:
            return model.get_grid_z(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_grid_connectivity(self, long_var_name):
        try:
            return model.get_grid_connectivity(long_var_name)
        except Exception as e:
            raise ModelException(str(e))
    
    def get_grid_offset(self, long_var_name):
        try:
            return model.get_grid_offset(long_var_name)
        except Exception as e:
            raise ModelException(str(e))

    # extended bmi functions

    def initialize_config(self, config_file):
        if config_file is None or config_file == "":
            config_file = None
        
        try:
            model.initialize_config(config_file)
        except Exception as e:
            raise ModelException(str(e)) 

    def initialize_model(self):
        try:
            model.initialize_model()
        except Exception as e:
            raise ModelException(str(e)) 

    def set_start_time(self, start_time):
        try:
            model.set_start_time(start_time)
        except Exception as e:
            raise ModelException(str(e)) 

    def set_end_time(self, end_time):
        try:
            model.set_end_time(end_time)
        except Exception as e:
            raise ModelException(str(e)) 

    def get_attribute_names(self):
        try:
            return model.get_attribute_names()
        except Exception as e:
            raise ModelException(str(e)) 

    def get_attribute_value(self, attribute_name):
        try:
            return model.get_attribute_value(attribute_name)
        except Exception as e:
            raise ModelException(str(e)) 

    def set_attribute_value(self, attribute_name, attribute_value):
        try:
            model.set_attribute_value(attribute_name, attribute_value)
        except Exception as e:
            raise ModelException(str(e)) 

    def save_state(self, destination_directory):
        try:
            model.save_state(destination_directory)
        except Exception as e:
            raise ModelException(str(e)) 

    def load_state(self, source_directory):
        try:
            model.load_state(source_directory)
        except Exception as e:
            raise ModelException(str(e)) 

    
def handleSIGINT(sig, frame):
    #clean up state or what ever is necessary
    sys.exit(0)

if __name__ == '__main__':

    #setup root logger to catch any warnings and errors. Most models will override and/or add to these settings.
    logging.basicConfig(level=logging.WARNING)
    
    model_module_name = sys.argv[1]
    model_class_name = sys.argv[2]
    
    model_module = importlib.import_module(model_module_name)
    model_class = getattr(model_module, model_class_name)
    
    model = model_class()
    
    handler = ModelHandler(model)
    processor = Processor(handler)
    
    transport = TSocket.TServerSocket(host=sys.argv[3],port=sys.argv[4])

    tfactory = TTransport.TBufferedTransportFactory()
    pfactory = TBinaryProtocol.TBinaryProtocolFactory()
    
    server = TServer.TSimpleServer(processor, transport, tfactory, pfactory)
    
    signal.signal(signal.SIGINT, handleSIGINT)
    
    logger.info(server)

    logger.info("serving")
    
    server.serve()
    
    logger.info("done")
