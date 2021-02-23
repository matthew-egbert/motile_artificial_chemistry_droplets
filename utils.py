import logging as logger
import color_logger
import sys

logging_parameters = {
    'level' : logger.INFO,
#    'format' : '%(levelname)s:%(filename)s:%(message)s'
}
logger.basicConfig(**logging_parameters) # filename='example.log',

#logger.setLevel(logger.DEBUG)

### Local Variables: ###
### python-main-file: "main.py" ###
### End: ###
