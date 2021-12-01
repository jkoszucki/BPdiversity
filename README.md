Determination of phage variants. 

1. Megablast all phages vs all phages.
Sometimes mutiple hits are given for two phages.

2. Hits prefiltering. Leave only significant hits.
pident >= 0.75
eval <= 10**-3
bitscore >= 100

3. Hits are merged and query and subject coverage is calculated. 

4. Phages are considered as one phage variant if (confi.yaml):
query coverage >= 0.9
subject coverage >= 0.9
pident >= 0.9
