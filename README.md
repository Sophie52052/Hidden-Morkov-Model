# (Hidden) Markov Model - Genome
###### tags: `Homework`

* **Data**
    * GRCh38, Fasta
    * NC_000006.12 Homo sapiens chromosome 6, GRCh38.p13
    * Bases [100,000–199,999]
    * Starting with ttggtaccat and ending in CTTTGCCTG.
    * Extract bases [100,000–199,999]; starting with ttggtaccat and ending in CTTTGCCTG.
* **Markov Chain**
    1. **Order 0**
        * 意義：
            * 考慮每個鹼基獨立出現在序列中的機率
        * 做法：
            * 統計A、T、C、G個別在此序列中出現的機率
            * 依據序列出現鹼基的順序把個別出現的機率成起來
            * 答案即為所求
        * 結果：
            * order 0 : **-197688.94757656846**

        * 例子：
             `A T C A T C G G A G`
             `A:0.3(P(A)) T:0.2(P(T)) C:0.2(P(C)) G:0.3(P(G))`
             `0.3 * 0.2 * 0.2 ....`
             
        
    
        

     
    2. Order 1
        * 意義：
            * 考慮每個鹼基以及其前一個鹼基出現在序列中的機率
        * 做法：
            * 統計AA、AT、AC、AG、TA、TC...(共4*4=16種)出現在此序列中的機率
            * 依據序列出現鹼基的順序把個別出現的機率成起來
            * 答案即為所求
        * 結果：
            * order 1 : **-193211.25299986434**
        * 例子：
             `A T C A T C G G A G`
             `AA:P(AA) AT:P(AT) AC:P(AC) AG:P(AG)...`
             `P(A) * P(AT) * P(TC) * P(CA)...`

    4. Order 2
        * 意義：
            * 考慮每個鹼基以及其前二個鹼基出現在序列中的機率
        * 做法：
            * 統計AAA、AAT、AAC、AAG、ATA、ATT、ATC...(共4* 4 *4=64種)出現在此序列中的機率
            * 依據序列出現鹼基的順序把個別出現的機率成起來
            * 答案即為所求
        * 結果：
            * order 2 : **-192098.4269535206**
        * 例子：
             `A T C A T C G G A G`
             `AAA:P(AA) AAT:P(AT) AAC:P(AC) AAG:P(AG)...`
             `P(A) * P(AT) * P(ATC) * P(TCA)...`

* **Hidde Markov Model**
    1. 使用**Expectation-Maximization Algorithm**：
        * 參考演算法筆記的演算法並作改良
        * http://www.csie.ntnu.edu.tw/~u91029/index.html
    2. 經歷：
        * 起初使用**Viterbi Algorithm**-但發現python會一直出現**遞迴太深**的錯誤
        * 於是改用**Expectation-Maximization Algorithm**-但也同時發現forward在計算時會算到-inf
        * 因此從**Expectation-Maximization Algorithm改良**
            * 由於這兩種演算法的問題就是在暴力解開並找最大值時，會遞太深
            * 因此，目前改良過的forward就不計算所有可能性，只在當下最大的機率下在去往下找
    3. Model建構：
        * 4個state
        * 每個state裡面有4個機率
        ![](https://i.imgur.com/pdVuN5x.png)
        * **轉移機率**設置為Markov Chain統計出的Order 2的機率
        * State內部機率為Markov Chain統計出的Order 0的機率
    4. 結果：100000-1999999 Maximun Probability :  **-346983.20505182043**
    5. 其他：
        * 找出額外100KB的資料-使用chromosome 6 [300,000–401,999]
        * 結果：300000-401999 Maximun Probability :  **-358268.38464508875**

* 運行結果圖

    ![](https://i.imgur.com/JKOOrhk.png)
    ![](https://i.imgur.com/zzQJRyw.png)



