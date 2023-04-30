INSTRUKCJA OBSŁUGI

1. DO CZEGO SŁUŻY PROGRAM?
Program służy do przeliczania współrzędnych. Obsługuje on następujące transformacje:
przeliczenie współrzędnych geocentrycznych XYZ na współrzędne geodezyjne BLH, 
przeliczenie współrzędnych geodezyjnych BLHdo współrzędnych prostokątnych w układzie 2000, 
przeliczenie współrzędnych geodezyjnych BLH do współrzędnych prostokątnych w układzie 1992, 
przeliczenie współrzędnych geodezyjnych BLH do współrzędnych geocentrycznych XYZ, 
przeliczenie współrzędnych geocentrycznych XYZ do współrzędnych topocentrcznych NEU, 
przeliczenie współrzędnych geodezyjnych BLH do współrzędnych geocentrycznych XYZ, 

Program obsługuje parametry modeli elipsoid: GRS80, WGS84, Krassowskiego

2. WYMAGANIA:
Aby program zadziałał na danym komputerze/laptopie musi on mieć zainstalowanego Pythona wersję 3.10, bibliotekę Numpy, bibliotekę Argparse, a także zintegrowane środowisko-Spyder wersja 5.1.5.

3. SYSTEM OPERACUJNY NA KOMPUTERZE A PROGRAM:
Program został napisany dla systemu operacyjnego Windows 10

4. UŻYTKOWANIE PROGRAMU:
Aby skorzystać z programu użytkownik powinien posiadać przygotowany wcześniej plik (rozszerzenie .txt) z danymi wyjściowymi zapisanymi w określony poniżej sposób:

a) W przypadku współrzędnych geocentrycznych XYZ jako danych wejściowych plik tekstowy powinien zawierać w kolejnych liniach owe 
współrzędne, zapisane w następującej kolejności: X, Y, Z, przy czym współrzędne musza być one oddzielone od siebie przecinkami, 
a części dziesiętne oddzielone od jedności kropką. 
Przykładowo: 1000.000, 1200.000, 1300.000
Dane wyjściowe zostaną zapisane jako: B(decimal), B(minutes), B(seconds), L(decimal), L(minutes), L(seconds), H (metry).

b) W przypadku współrzędnych geodezyjnych jako danych wejściowych plik tekstowy powiniem zawierać w kolejnych liniach owe współrzędne 
zapisane w następującej kolejności: B(decimal), B(minutes), B(seconds), L(decimal), L(minutes), L(seconds), H przy czym współrzędne musza być one oddzielone od siebie przecinkami, a części dziesiętne oddzielone od jedności kropką.
Przykładowo: 52,16,22.000,17,11,21.888,131.000
W przypadku transformacji do współrzędnych geocentrycznych dane wyjściowe zostaną zapisane jako współrzędne geocentryczne: X,Y,Z, 
przy czym każda ze współrzędnych określona do części tysięcznych. 
W przypadku transformacji do układu 1992 dane wyjściowe zostaną zapisane jako współrzędne geocentryczne: X1992,Y1992 
przy czym każda ze współrzędnych określona do części tysięcznych.
W przypadku transformacji do układu 2000 dane wyjściowe zostaną zapisane jako współrzędne geocentryczne: X2000,Y2000 
przy czym każda ze współrzędnych określona do części tysięcznych.


c) W przypadku współrzędnych geocentrycznych XYZ początka odcinka i współrzędnych geocentrycznych XYZ końca odcinka jako danych wejściowych
plik tekstowy powinien zawierać w kolejnych liniach owe współrzędne zapisane w następującej kolejności:
Xpoczątkowe, Ypoczątkowe, Zpoczątkowe, Xkońcowe, Ykońcowe, Zkońcowe  przy czym współrzędne musza być one oddzielone od siebie przecinkami,
a części dziesiętne oddzielone od jedności kropką. 
Przykładowo:  1000.000, 1200.000, 1300.000, 1001.000, 1202.000, 1303.000

Aby skorzystać z programu i danych zapisanych na pliku wejściowym, wywołujemy go przez wiersz poleceń. Program korzysta z biblioteki argparse, 
zatem przy wowołaniu podajemy argumenty. Argument "-m" oznacza model elipsoidy, na której dokonujemy przeliczeń współrzędnych, przy czym program obsługuje trzy, które użytkownik ma do wyboru: "wgs84","grs80", "krassowski".
"-t" oznacza transformację, z jakiej użytkownik chce skorzystać, przy czym ma do wyboru: "Hirvonen", "BLHtoXYZ", "BLHto2000", "BLHto1992".
Argument ,,-from_file" oznacza plik z którego program pobiera dane wejściowe dla danej transformacji, a argument "-to_file" plik do którego zostaną zapisane dane wyjściowe. Przy owych polach jako argument użytkownik musi podać ścieżkę/nazwę pliku.

