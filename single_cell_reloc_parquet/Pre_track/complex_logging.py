#%%
import logging
import os
from joblib import Parallel, delayed
from multiprocessing import Queue, current_process
from logging.handlers import QueueHandler, QueueListener
from time import sleep

def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    logging.getLogger().setLevel(logging.INFO) #* Logger 1

    pyspark_log = logging.getLogger('pyspark') #*Logger 2
    pyspark_log.setLevel(logging.WARNING)

    logging.getLogger("py4j").setLevel(logging.ERROR) #*Logger 3
    logging.getLogger("py4j.java_gateway").setLevel(logging.ERROR)

    # Create a queue handler and listener for logging in child processes
    log_queue = Queue(-1)
    queue_handler = QueueHandler(log_queue)
    queue_listener = QueueListener(log_queue, logging.getLogger())

    # Start the queue listener in a separate thread
    queue_listener.start()

def calculate_square(number, log_queue):
    sleep(1)  # Simulate a time-consuming task
    result = number ** 2
    logging.info(f"[{current_process().name}] Square of {number} is {result}")
    log_queue.put(f"[{current_process().name}] Square of {number} is {result}")
    return result

setup_logging()
logging.info(f"Starting multiprocessing")
# List of numbers to calculate squares
numbers = [1, 2, 3, 4, 5]

# Using joblib for multiprocessing
n_jobs = os.cpu_count()

# Create a queue for logging in child processes
log_queue = Queue(-1)

# Create a delayed function with the log_queue as an argument
def delayed_func(num):
    return delayed(calculate_square)(num, log_queue)

results = Parallel(n_jobs=n_jobs)(delayed_func(num) for num in numbers)

# Stop the queue listener once all child processes have finished logging
queue_listener.stop()

print(f"Squares: {results}")
# %%
