# Spustenie:
Python - `python test.py`(windows) <br>
C - `comp.bat`(windows) <br>
# Progress:

## 23.4.2023:
Upravil som funckie encode/decode, podľa [draftu ristretto255](https://datatracker.ietf.org/doc/draft-irtf-cfrg-ristretto255-decaf448/).  <br>
Taktiež som upravil funkciu feq() pretože na vracala chybne výsledky.  <br>
Prečistil som kód od warningov (ostali len tie s unused variables - tie casto pouzivam pri debuggovani takze som ich tam este nechal). <br>  

<br>

Zaver:
Dokončil som Impl. ristretta255 pre Python ale aj C, výsledky dávaju identické a ďalším krokom je otestovanie, či impl. funguju správne. <br>


## 28.4.2023:
Upravil som bug v decode pre C aj Python kod <br>
Pridal som "testovaci vektor", t.j. validny vektor u8 VECTOR_TEST2[32] v maine.  <br>
Warningy som ponechal (ostali len tie s unused variables - tie casto pouzivam pri debuggovani takze som ich tam este nechal). <br>  

<br>

Zaver:
Otestoval som decoding/encoding na vektore VECTOR_TEST2. Po kombinácii de/en-code dostanem input==output. <br>
Z toho súdim, že impl. by mala byť zatiaľ v poriadku. Ďalším krokom je porovnanie s dalsimi impl. a vytvorenie hash_to_group (kt. chyba v drafte) <br>
