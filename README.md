# TransposableElements
Call TEs with PopoolationTE2 and remove low complexity regions. For PopTE2 make sure the TE library is good (RepBase TEs) and remove TEs that blast to genes or show strange coverage. Also make sure the TEs are >400bp <br/>
Step1: Filter to make bed file within an individuals <br/>
Step2: Merge bed file in a population<br/>
Step3: Check TE calls using the discordant short read pairs to remove PCR artifacts and low coverage calls <br/>
Step4: Check TE calls using the chimeric reads <br/>
