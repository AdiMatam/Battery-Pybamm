M = 3
N = 10

counts = {}
for i in range(M):
    for j in range(N):
        counts[(i, j)] = 0


## IAPP check

## parallel nodes
for j in range(N-1):
    counts[(M-1, j)] += 1
counts[(M-1,N-1)] += 1

## series nodes
for n in range(N):
    for m in range(M-1):
        counts[(m, n)] += 1


## PHI CHECK

# parallel
for n in range(N-1):
    counts[(0, n)] += 1
    counts[(M-1,n)] += 1
counts[(0, N-1)] += 1
counts[(M-1,N-1)] += 1

# series
for n in range(N):
    for m in range(M-1):
        counts[(m, n)] += 1
        counts[(m+1, n)] += 1


## show counts
for key, value in counts.items():
    print(f"{key} -> {value}")
