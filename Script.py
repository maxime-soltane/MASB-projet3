from pybloom_live import ScalableBloomFilter   
    
    
def init_bloom_filter(sequence, k) -> ScalableBloomFilter:
    """ 
    Initializes a Bloom filter with all k-mers in the genome

    Returns:
    The Bloom filter containing the genome's kmers
    """
    #Initialize the bloom filter with an estimated capacity and the error rate of the filter returning false positives
    bloom = ScalableBloomFilter(initial_capacity=len(sequence)//k * 2, error_rate=0.01)
        
    # Loop through the sequence to extract and add kmers in the bloom filter
    for i in range(len(sequence) - k + 1):
        bloom.add(sequence[i:i+k])

    return bloom

