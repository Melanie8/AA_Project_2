import time
import datetime as dt

print dt.datetime.today()
print dt.datetime.today().timestamp()
start_time = dt.datetime.today().timestamp()
i = 0
while(True):
    time.sleep(0.1)
    time_diff = dt.datetime.today().timestamp() - start_time
    i += 1
    print(i / time_diff)
