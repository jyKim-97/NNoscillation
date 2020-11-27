from functools import partial
import multiprocessing
import time


def singleCount(name):
    cnt = 0
    for i in range(int(1e6)):
        cnt += 1
        if (i%int(1e5) == 0):
            print(name, ":", cnt)
        
# names = ['1', '2', '3', '4']

names = ['%d'%(i) for i in range(100)]

tic = time.time()
pool = multiprocessing.Pool(processes=12)

# cnt = 0

# func = partial(singleCount, cnt)
pool.map(singleCount, names)
pool.close()
pool.join()
toc = time.time()

tic_single = time.time()
for name in names:
    singleCount(name)

toc_single = time.time()

print('MULTI time = %.3fs\n'%(toc-tic))
print('MULTI time = %.3fs\n'%(toc_single-tic_single))