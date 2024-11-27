#%%
import logging
from logging.handlers import SysLogHandler
import time
from sys import _getframe
from datetime import datetime

# def setup_logging(logger_name):
# 	global PAPERTRAIL_HOST
# 	global PAPERTRAIL_PORT
# 	global formatter
# 	global file_handler
# 	global send_handler

# 	exec(f"global {logger_name}")
# 	PAPERTRAIL_HOST = "logs.papertrailapp.com"
# 	PAPERTRAIL_PORT = 31634

# 	formatter = logging.Formatter("%(asctime)s %(funcName)s %(levelname)s %(message)s %(lineno)d")
# 	logging.Formatter.converter = time.gmtime #* Set the date to be gmtime instead of local

# 	file_handler =logging.FileHandler(f'{logger_name}.log')
# 	file_handler.setFormatter(formatter)
# 	file_handler.setLevel(logging.warning)

# 	send_handler = SysLogHandler(address=(PAPERTRAIL_HOST, PAPERTRAIL_PORT))
# 	send_handler.setFormatter(formatter)
# 	send_handler.setLevel(logging.ERROR)

# 	exec(f"{logger_name} = logging.getLogger('{logger_name}')")
# 	exec(f"{logger_name}.addHandler(send_handler)")
# 	exec(f"{logger_name}.addHandler(file_handler)")
# 	return(None)

def setup_logging(logger_name, microfluidics_results):
	PAPERTRAIL_HOST = "logs.papertrailapp.com"
	PAPERTRAIL_PORT = Secret

	log_folder = os.path.join(microfluidics_results, 'logs')
	try:
		os.mkdir(log_folder)
	except FileExistsError:
		pass

	formatter = logging.Formatter("%(asctime)s %(funcName)s %(levelname)s %(message)s %(lineno)d")
	logging.Formatter.converter = time.gmtime #* Set the date to be gmtime instead of local
	now = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
	file_handler =logging.FileHandler(f'{log_folder}/{logger_name}_{now}.log')
	file_handler.setFormatter(formatter)
	file_handler.setLevel(logging.WARNING)

	send_handler = SysLogHandler(address=(PAPERTRAIL_HOST, PAPERTRAIL_PORT))
	send_handler.setFormatter(formatter)
	send_handler.setLevel(logging.ERROR)

	logger = logging.getLogger(f'{logger_name}')
	logger.addHandler(send_handler)
	logger.addHandler(file_handler)
	return(logger)


#%%

# # def define_logger():
# logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)
# logger.fo



# logger.addHandler(file_handler)
# # logger.addHandler(send_handler)

# logger.critical('Critical Error inside')
# 	# return(logger)

# # logger = logger_setup()

# def test_function():
# 	logger = logging.getLogger(__name__)
# 	logger.critical("New_message")



# def main() -> None:
#     logger = logging.getLogger("arjancodes")
#     logger.setLevel(logging.DEBUG)
#     handler_send = SysLogHandler(address=(PAPERTRAIL_HOST, PAPERTRAIL_PORT))
#     logger.addHandler(handler_send)

#     logger.warning("This is a debug message.")


# def main() -> None:
# 	hostname=socket.gethostname()
# 	IPAddr=socket.gethostbyname(hostname)
# 	print("Your Computer Name is:"+hostname)
# 	print("Your Computer IP Address is:"+IPAddr)
# 	logging.basicConfig(
# 		level=logging.DEBUG,
# 		format="%s %(asctime)s %(levelname)s %(message)s, %{hostname}",
# 		datefmt="%Y-%m-%d %H:%M:%S",
# 	)

# 	logging.debug("This is a debug message.")
# 	logging.info("This is an info message.")
# 	logging.warning("This is a warning message.")
# 	logging.error("This is an error message.")
# 	logging.critical("This is a critical message.")


# if __name__ == "__main__":
# 	main()
# %%
