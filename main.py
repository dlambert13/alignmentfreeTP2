from loading import load_directory
from kmers import stream_kmers, kmer2str
from collections import Counter
import timeit
import heapq
import time


def index_file(sequences, k):
    return Counter(list_kmers(sequences, k))


def list_file(sequences, k):
    lst = []

    # Construct the list of all the kmers of all the sequences
    for seq in sequences:
        lst.extend([min(x) for x in stream_kmers(seq, k)])

    # Sort the list
    lst.sort()

    return lst


def partitionminhash_sketch():
    pass


def intersect_index(index, sequences, k):
    """ Create an index containing all the kmers of the input sequences
    """
    index_uniq = index.copy()
    query_uniq = 0
    intersection = 0

    for seq in sequences:
        for kmer, rkmer in stream_kmers(seq, k):
            minmer = min(kmer, rkmer)

            # Query not in index
            if minmer not in index_uniq:
                query_uniq += 1
            # Query in index => intersection
            else:
                intersection += 1
                index_uniq[minmer] -= 1
                if index_uniq[minmer] == 0:
                    del index_uniq[minmer]

    return sum(index_uniq.values()), intersection, query_uniq


def intersect_sorted_lists(lst1, lst2):
    idx1 = idx2 = 0
    # Dataset specific or intersection kmer counts
    A = inter = B = 0

    while (idx1 < len(lst1) and idx2 < len(lst2)):
        kmer1 = lst1[idx1]
        kmer2 = lst2[idx2]

        # Same kmer => intersection
        if kmer1 == kmer2:
            inter += 1
            idx1 += 1
            idx2 += 1
        # first list specific
        elif kmer1 < kmer2:
            A += 1
            idx1 += 1
        # second list specific
        else:
            B += 1
            idx2 += 1

    # Add remaining kmers of the non empty list
    A += len(lst1) - idx1
    B += len(lst2) - idx2

    return A, inter, B


def list_kmers(sequences, k):
    kmers = []
    for seq in sequences:
        kmers.extend([min(kmer, rkmer) for kmer, rkmer in stream_kmers(seq, k)])
    return kmers


def similarity(A, inter, B):
    # +1 added for pseudocount. Avoid divisions by 0
    A_similarity = inter / (inter + A + 1)
    B_similarity = inter / (inter + B + 1)

    return A_similarity, B_similarity


def jaccard(A, inter, B):
    return inter / (A + inter + B)


def maxheap(kmer_list):
    list_copy = [-kmer for kmer in kmer_list]
    heapq.heapify(list_copy)
    return list_copy


def lightmax(input_list):
    buffer_list = maxheap(input_list)
    return - buffer_list[0]


def xorshift64(i):
    # after https://en.wikipedia.org/wiki/Xorshift
    i ^= i << 13
    i ^= i >> 7
    i ^= i << 17
    return i


def minhash_sketch(sequences, k, s):
    l = []
    for seq in sequences:
        for km, rkm in stream_kmers(seq, k):
            kmer = min(xorshift64(km), xorshift64(rkm))

            if len(l) < s:
                l.append(kmer)
                if len(l) == s: heapq.heapify(l)

            else:
                list_max = lightmax(l)
                if kmer < list_max:
                    heapq.heappop(l)
                    heapq.heappush(l, kmer)
    
    l.sort()
    return l


if __name__ == "__main__":
    """
    k = 3
    s = 5
    test_seqs_a = ["ACCACG", "TACCAC", "TCCACT", "TTTACA", "ACCATG"]
    test_seqs_b = ["ACCTCG", "TATTAC", "TCCACT", "TTAACA", "ACCATG"]
    mhs_a = minhash_sketch(test_seqs_a, s, k)
    mhs_b = minhash_sketch(test_seqs_b, s, k)
    print(mhs_a, mhs_b)
    A, inter, B = intersect_sorted_lists(mhs_a, mhs_b)
    print(similarity(A, inter, B), " ; ", jaccard(A, inter, B))
    """
    # original pre-modification "tp2" script
    # Load all the files in a dictionary
    files = load_directory("data")

    k = 21
    s = 1000
    # Loading
    indexes = {f:index_file(files[f], k) for f in files}
    lists = {f:list_file(files[f], k) for f in files}
    minhash_sketches = {f:minhash_sketch(files[f], k, s) for f in files}

    filenames = list(files.keys())
    #for i in range(len(files)):
    for i in [0]:
        for j in [1]:
        #for j in range(i+1, len(files)):
            # Method 1 using index
            print("index")
            A, inter, B = intersect_index(indexes[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
            # Method 2 using sorted lists
            print("lists")
            A, inter, B = intersect_sorted_lists(lists[filenames[i]], lists[filenames[j]])
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
            # Method 3 using minhash sketches
            print("minhash")
            A, inter, B = intersect_sorted_lists(
                minhash_sketches[filenames[i]],
                minhash_sketches[filenames[j]]
                )
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
