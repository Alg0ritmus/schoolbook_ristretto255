# Spustenie:
Python - `python test.py`(windows) <br>
C - `comp.bat`(windows) <br>
# Progress:

## 14.6.2023
Doplnil som automatizovane jednotkove testy v maine a taktiez som pridal prepinac na vypis jednotlivych testov (-DVERBOSE_FLAG) spolu <br>
s prepinacom na chybove hlasky, ktore sa vyskytuju v pripade roznych chyb ktore mozu nastat napr. pri ristretto255_decode (-DDEBUG_FLAG).<br>
Prepinace (v subore comp.bat) by mohli vyzerat nejak takto: <br>
1) `gcc -O0 -DDEBUG_FLAG  -Wall -Wno-unused-value -Wextra -o main main.c ristretto255.c && main.exe`<br>
2)`gcc -O0 -DDEBUG_FLAG  -Wall -Wno-unused-value -Wextra -o main main.c ristretto255.c && main.exe`<br>
3)`gcc -O0 -DVERBOUS_FLAG -DDEBUG_FLAG  -Wall -Wno-unused-value -Wextra -o main main.c ristretto255.c && main.exe`<br>

Co sa tyka CT(constant time) funkcii:<br>
1) upravil som funkciu fabsolute, tak aby pracovala v konstantnom case<br>
2) CT_SELECT je vlastne swap25519 ked sa nad tym zamyslime, cize vyuzil som funkciu swap25519 priamo z tweetnacl<br>

 <br>


## 10.6.2023
Uspesne som implementoval som funkciu hash_to_group a otestoval pomocou testovacich vektorov z draftu aj pre jazyk C. Trochu som kod precistil <br>
Dalsim krokom bude re-implementovat niektore funkcie ako vykonavali ukony v konstantnom case (napr absolut. hodnota + implementovat (CT_SELECT, CT_EQ, CT_NEG, CT_ABS)).<br>

 <br>

## 24.5.2023
Uspesne som implementoval som funkciu hash_to_group a otestoval pomocou testovacich vektorov z draftu. Od poslednej upravy kodu som prisiel na to, <br>
ze chyba v hash_to_group suvisela najma s maskovanim LSB + s nejakymi inymi mensimi chybami v kode. Zaroven scitanie bodov bodov vo funkcii hash_to_group <br>
vyuziva scinatie edwardsovych bodov z libsodia (**libsodium_add.py**). Dalej som precistil a okomentoval Python kod a pridal testovacie vektory (z draftu/ rust impl.). <br>
Dalsim krokom bude implementovat hash_to_group aj v jazyku C a zaroven implementovat aj nasobenie bodu ristretto skalarom.  <br>

 <br>

## 13.5.2023
Otestoval som ci kniznica [ge25519](https://pypi.org/project/ge25519/) dava vysledky zhodujuce sa s testovacimi <br>
vektormi v drafte (ano, zhoduju sa). Kniznica je vhodna na detailnejsiu analyzu v pripade potreby.<br>
Bohuzial kvoli inej reprezentacii cisel (oni pouzivaju [fe25519](https://pypi.org/project/fe25519/)) neviem otestovat vysledky v medzivypoctoch (iba konecny vysledok). <br>
Funkciu `hash_to_group()` sa mi zatial nepodarilo uspesne implementovať, pravdepodobne bude problem v scitani.<br>
Prave kvoli tomu som vyskusal viacero sposobov scitania bodu na Edwardsovej krivke (vid. nizsie), avsak zatial neuspesne.<br>
Precistil som kod a taktiež ho rozdelil kod do viacerych suborov kvoli prehladnosti:<br>
  1) constants.py - obsahuje vsetky konstanty (z draftu) <br>
  2) scitanie bodu: <br>
      2.1) ADD z TweetNaCl sa nachádza v subore: **bernstein.py**, názov funkcie je **Fadd()** <br>
      2.2) scitanie z libsodia sa nachádza v subore: **libsodium_add.py**, názov funkcie je **ge25519_p3_add()** <br>
      2.3) scitanie z [paper](https://eprint.iacr.org/2008/522.pdf) sa nachádza v subore: **edwardsc_paper_add.py**, názov funkcie je **edwards_addition()** <br>
 <br>

## 28.4.2023:
Upravil som bug v decode pre C aj Python kod <br>
Pridal som "testovaci vektor", t.j. validny vektor u8 VECTOR_TEST2[32] v maine.  <br>
Warningy som ponechal (ostali len tie s unused variables - tie casto pouzivam pri debuggovani takze som ich tam este nechal). <br>  

<br>

Zaver:
Otestoval som decoding/encoding na vektore VECTOR_TEST2. Po kombinácii de/en-code dostanem input==output. <br>
Z toho súdim, že impl. by mala byť zatiaľ v poriadku. Ďalším krokom je porovnanie s dalsimi impl. a vytvorenie hash_to_group (<s>kt. chyba v drafte?</s> <i>element derivation function</i>) <br>

## 23.4.2023:
Upravil som funckie encode/decode, podľa [draftu ristretto255](https://datatracker.ietf.org/doc/draft-irtf-cfrg-ristretto255-decaf448/).  <br>
Taktiež som upravil funkciu feq() pretože na vracala chybne výsledky.  <br>
Prečistil som kód od warningov (ostali len tie s unused variables - tie casto pouzivam pri debuggovani takze som ich tam este nechal). <br>  

<br>

Zaver:
Dokončil som Impl. ristretta255 pre Python ale aj C, výsledky dávaju identické a ďalším krokom je otestovanie, či impl. funguju správne. <br>




