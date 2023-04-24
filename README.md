# Progress:

## 23.4.2023:
Upravil som funckie encode/decode, podľa [draftu ristretto255](https://datatracker.ietf.org/doc/draft-irtf-cfrg-ristretto255-decaf448/).  <br>
Taktiež som upravil funkciu feq() pretože na vracala chybne výsledky.  <br>
Prečistil som kód od warningov (ostali len tie s unused variables - tie casto pouzivam pri debuggovani takze som ich tam este nechal). <br>  

<br>

Zaver:
Dokončil som Impl. ristretta255 pre Python ale aj C, výsledky dávaju identické a ďalším krokom je otestovanie, či impl. funguju správne. <br>
