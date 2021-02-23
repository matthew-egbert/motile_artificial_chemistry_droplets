import os
from datetime import datetime
import time
import shutil
import jinja2

class EventLogger(object):
    def __init__(self):
        self.log_path = 'log'    
        shutil.rmtree(self.log_path)
        os.mkdir(self.log_path)

        self.jinja_loader = jinja2.FileSystemLoader(searchpath="./log_templates")
        self.jinja_env = jinja2.Environment(loader=self.jinja_loader)
        
        events = []

        self.avalanche_count = 0

    def record_avalanche(self):
        shutil.copy('graph.png',
                    os.path.join(self.log_path,f'graph_post_avalanche_{self.avalanche_count}.png'))
        self.avalanche_count += 1
        pass

    def close(self):
        template = self.jinja_env.get_template('index.html')
        out = template.render()  # this is where to put args to the template renderer
        print(out)
        

    
        
