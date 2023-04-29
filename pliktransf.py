import numpy as np

class Transformacje:
    
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "krassowski":
            self.a = 6378245.00000
            self.b = 6356863.01877
        else:
            raise NotImplementedError(f"{model} niezaimplementowany model")
        self.flat = (self.a - self.b) / self.a
        self.ecc = np.sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
        
        
    def Hirvonen(self, X, Y, Z):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne geocentryczne, 

        Returns
        -------
        phi
            [stopnie dziesiętne] - szerokość geodezyjna
        lam
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        p   = np.sqrt(X**2 + Y**2)           # promień
        phi_poprzednie= np.arctan(Z / (p * (1 - self.ecc2)))    
        phi = 0
        while abs(phi_poprzednie - phi) > 0.000001/206265:    
            phi_poprzednie = phi
            N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi_poprzednie)**2)
            h = p / np.cos(phi_poprzednie) - N
            phi = np.arctan((Z/p) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lam = np.arctan(Y/X)
        N = self.a / np.sqrt(1 - self.ecc2 * (np.sin(phi))**2);
        h = p / np.cos(phi) - N       
        
        return np.rad2deg(phi), np.rad2deg(lam), h 
         
            
    def BLHto2000(self, phi, lam, h):
            '''
            
            Transformacja współrzędnych geodezyjnych φ, λ, h
            na układ współrzędnych prostokątnych płaskich w układzie 2000 
            to przekształcenie z jednego układu odniesienia na drugi. 
            Wartości φ i λ wyrażają położenie punktu na powierzchni elipsoidy referencyjnej, a układ 2000 jest planimetrycznym układem współrzędnych prostokątnych płaskich,
            który stosuje się do map, a jego wartości współrzędnych wyrażone są w metrach.
            Parameters
            ----------
            phi, lam, h : FLOAT
                 współrzędne geodezyjne, 

            Returns
            -------
            X2000
                [metry] - współrzędna prostokątna płaska X
            Y2000
                [metry] - współrzędna prostokątna płaska Y 
           
            '''
            
            if lam>13.5 and lam <16.5:
                s = 5
                l0 = 15.00000
            elif lam>16.5 and lam <19.5:
                s = 6
                l0 = 18.00000
            elif lam>19.5 and lam <22.5:
                s = 7
                l0 = 21.00000
            elif lam>22.5 and lam <25.5:
                s = 8
                l0 = 24.00000
            
            l0=np.deg2rad(l0)
            phi=np.deg2rad(phi)
            lam=np.deg2rad(lam)
            
            a2 = self.a**2
            b2 = a2 * (1 - self.ecc2)
            e_2 = (a2 - b2)/b2
            dl = lam - l0
            dl2 = dl**2
            dl4 = dl**4
            t = np.tan(phi)
            t2 = t**2
            t4 = t**4
            n2 = e_2 * (np.cos(phi)**2)
            n4 = n2 ** 2
            N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
            e4 = self.ecc2**2
            e6 = self.ecc2**3
            A0 = 1 - (self.ecc2/4) - ((3*e4)/64) - ((5*e6)/256)
            A2 = (3/8) * (self.ecc2 + e4/4 + (15*e6)/128)
            A4 = (15/256) * (e4 + (3*e6)/4)
            A6 = (35*e6)/3072
            sigma = self.a * ((A0 * phi) - A2 * np.sin(2*phi) + A4 * np.sin(4*phi) - A6 * np.sin(6*phi))
            xgk = sigma + ((dl**2)/2) * N * np.sin(phi) * np.cos(phi) * (1 + ((dl**2)/12)*(np.cos(phi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(phi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
            ygk = dl * N * np.cos(phi) * (1 + (dl2/6) * (np.cos(phi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(phi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
            
            X2000 = xgk * 0.999923
            Y2000 = ygk * 0.999923 + s * 1000000 + 500000
            return X2000, Y2000
        
        
        
    def BLHto1992(self, phi, lam, h):
            '''
        
            Transformacja współrzędnych geodezyjnych φ, λ, h
            na układ współrzędnych prostokątnych w układzie 1992 
            Wartości φ i λ wyrażają położenie punktu na powierzchni elipsoidy referencyjnej, a układ 1992 jest planimetrycznym układem współrzędnych prostokątnych płaskich,
            który stosuje się do map, a jego wartości współrzędnych wyrażone są w metrach.
            PARAMETERS:
            -------
            phi, lam, h : FLOAT
                 współrzędne geodezyjne, 

            Returns
            -------
            X2000
                [metry] - współrzędna prostokątna płaska X
            Y2000
                [metry] - współrzędna prostokątna płaska Y 
            '''
            phi=np.deg2rad(phi)
            lam=np.deg2rad(lam)
            #południk osiowy dla Polski równy:
            l0=19.00000
            l0=np.deg2rad(l0)
        
            a2 = self.a**2
            b2 = a2 * (1 - self.ecc2)
            e_2 = (a2 - b2)/b2
            dl = lam - l0
            dl2 = dl**2
            dl4 = dl**4
            t = np.tan(phi)
            t2 = t**2
            t4 = t**4
            n2 = e_2 * (np.cos(phi)**2)
            n4 = n2 ** 2
            N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
            e4 = self.ecc2**2
            e6 = self.ecc2**3
            A0 = 1 - (self.ecc2/4) - ((3*e4)/64) - ((5*e6)/256)
            A2 = (3/8) * (self.ecc2 + e4/4 + (15*e6)/128)
            A4 = (15/256) * (e4 + (3*e6)/4)
            A6 = (35*e6)/3072
            sigma = self.a * ((A0 * phi) - A2 * np.sin(2*phi) + A4 * np.sin(4*phi) - A6 * np.sin(6*phi))
            xgk = sigma + ((dl**2)/2) * N * np.sin(phi) * np.cos(phi) * (1 + ((dl**2)/12)*(np.cos(phi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(phi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
            ygk = dl * N * np.cos(phi) * (1 + (dl2/6) * (np.cos(phi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(phi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))

            X1992 = xgk * 0.9993 - 5300000
            Y1992 = ygk * 0.9993 + 500000
            return X1992, Y1992
        
    #def dXYZ2NEU(self, xp, yp, zp, xk, yk, zk):
        #'''
        #XYZ -> NEUp -
        #Transformacja z układu współrzędnych geocentrcznych XYZ 
        #na układ współrzędnych topograficznych NEU (układ lokalny związany z punktem odniesienia), 
        #służy do dokładnego okreslenia położenia punktu w lokalnym układzie współrzędnych.
        #PARAMETERS
        #-------
        #xp, yp, zp : FLOAT
        #     współrzędne geocentryczne punktu początkowego
        #xk, yk, zk : FLOAT
        #     współrzędne geocentryczne punktu końcowego     

        #Returns
        #-------
        #X2000
        #    [metry] - współrzędna prostokątna płaska X
        #Y2000
        #    [metry] - współrzędna prostokątna płaska Y 


        #'''
        #współrzędne początku  φ λ -> XYZ 
        #N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
        
        #obliczenie współrzędnych XYZ
        #dxyz = np.array([xk,yk,zk]) - np.array([xp,yp,zp])
        #X, Y, Z = dxyz - np.array([X0, Y0, Z0])
    
        # Macierz obrotu
        #sin_phi = np.sin(phi)
        #cos_phi = np.cos(phi)
        #sin_lam = np.sin(lam)
        #cos_lam = np.cos(lam)
    
        #R = np.array([[-sin_lam,           cos_lam,           0],
                      #[-sin_phi*cos_lam, -sin_phi*sin_lam,  cos_phi],
                      #[ cos_phi*cos_lam,  cos_phi*sin_lam,  sin_phi]])
    
        # Obliczenie współrzędnych NEU
        #NEU = np.dot(R, np.array([X, Y, Z]))
    
     
        #N = NEU[0]
        #E = NEU[1]
        #U = -NEU[2]
    
    
        #phi_stopnie = np.degrees(phi)
        #lam_stopnie = np.degrees(lam)
    
        #return (N, E, U, phi_stopnie, lam_stopnie)

if __name__=='__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='Transformacje współrzędnych')
    parser.add_argument('-m', type=str, help='Model elipsoidy', choices=['wgs84', 'grs80', 'krassowski'] )
    parser.add_argument('-t', type=str, help='Rodzaj transformacji', choices=['Hirvonen','BLH21992', 'BLHto2000', 'dXYZtoNEU'])
    parser.add_argument('-from_file', type=str, help='Ścieżka do pliku/nazwa pliku z którego pobieramy dane')
    parser.add_argument('-to_file', type=str, help='Ścieżka do pliku/nazwa pliku do którego zapisujemy dane')
    
    args = parser.parse_args()
    
    print(args.m, args.t, args.from_file, args.to_file)
    
    if args.m == 'grs80':
        if args.t == 'Hirvonen':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                X=[]
                Y=[]
                Z=[]
                for line in lines:
                    if "geocentryczne" in line:
                        continue
                    else:
                        rozdzielone_wsp=line.split(',')
                        X.append(float(rozdzielone_wsp[0]))
                        Y.append(float(rozdzielone_wsp[1]))
                        Z.append(float(rozdzielone_wsp[2]))
        
                for (x,y,z) in zip(X,Y,Z):
                    #utworzy obiekt o tych współrzędnych
                    geocentryczne_grs80=Transformacje(model = "grs80")
                    phi,lam,h = geocentryczne_grs80.Hirvonen(float(x), float(y), float(z))
                    wynik_transformacji.write(phi,lam,h)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
        elif args.t =='BLHto2000':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                B=[]
                L=[]
                H=[]
                for line in lines:
                    if "geodezyjne" in line:
                        continue
                    else:
                        rozdzielone_wsp = line.split(',')
                        Bd = rozdzielone_wsp[0]
                        Bm = rozdzielone_wsp[1]
                        Bs = rozdzielone_wsp[2]
                        
                        Bdeg=Bd+Bm/60+Bs/3600
                        B.append(Bdeg)
                        
                        Ld = rozdzielone_wsp[3]
                        Lm = rozdzielone_wsp[4]
                        Ls = rozdzielone_wsp[5]

                        Ldeg=Ld+Lm/60+Ls/3600
                        L.append(Ldeg)
                        
                        wys = rozdzielone_wsp[6]
                        H.append(wys)
                        
                for (b,l,h) in zip(B,L,H):
                    #utworzy obiekt o tych współrzędnych
                    geodezyjne2000_grs80=Transformacje(model = "grs80")
                    X2000, Y2000 = geodezyjne2000_grs80.BLHto2000(float(b), float(l), float(h))
                    print(f'Wynik transformacji BLH do układu 2000: X2000{X2000:.3f}, Y2000{Y2000:.3f}')
                    wynik_transformacji.write(X2000, Y2000)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
        elif args.t =='BLHto1992':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                B=[]
                L=[]
                H=[]
                for line in lines:
                    if "geodezyjne" in line:
                        continue
                    else:
                        rozdzielone_wsp = line.split(',')
                        Bd = rozdzielone_wsp[0]
                        Bm = rozdzielone_wsp[1]
                        Bs = rozdzielone_wsp[2]
                        
                        Bdeg=Bd+Bm/60+Bs/3600
                        B.append(Bdeg)
                        
                        Ld = rozdzielone_wsp[3]
                        Lm = rozdzielone_wsp[4]
                        Ls = rozdzielone_wsp[5]

                        Ldeg=Ld+Lm/60+Ls/3600
                        L.append(Ldeg)
                        
                        wys = rozdzielone_wsp[6]
                        H.append(wys)
                        
                for (b,l,h) in zip(B,L,H):
                    #utworzy obiekt o tych współrzędnych
                    geodezyjne1992_grs80=Transformacje(model = "grs80")
                    X1992, Y1992 = geodezyjne1992_grs80.BLHto1992(float(b), float(l), float(h))
                    print(f'Wynik transformacji BLH do układu 1992: X1992{X1992:.3f}, Y1992{Y1992:.3f}')
                    wynik_transformacji.write(X1992, Y1992)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
                
    elif args.m == 'wgs84':
        if args.t == 'Hirvonen':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                X=[]
                Y=[]
                Z=[]
                for line in lines:
                    if "geocentryczne" in line:
                        continue
                    else:
                        rozdzielone_wsp=line.split(',')
                        X.append(rozdzielone_wsp[0])
                        Y.append(rozdzielone_wsp[1])
                        Z.append(rozdzielone_wsp[2])
        
                for (x,y,z) in zip(X,Y,Z):
                    #utworzy obiekt o tych współrzędnych
                    geocentryczne_wgs84=Transformacje(model = "wgs84")
                    phi,lam,h = geocentryczne_wgs84.Hirvonen(float(x), float(y), float(z))
                    wynik_transformacji.write(phi,lam,h)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
        elif args.t =='BLHto2000':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                B=[]
                L=[]
                H=[]
                for line in lines:
                    if "geodezyjne" in line:
                        continue
                    else:
                        rozdzielone_wsp = line.split(',')
                        Bd = rozdzielone_wsp[0]
                        Bm = rozdzielone_wsp[1]
                        Bs = rozdzielone_wsp[2]
                        
                        Bdeg=Bd+Bm/60+Bs/3600
                        B.append(Bdeg)
                        
                        Ld = rozdzielone_wsp[3]
                        Lm = rozdzielone_wsp[4]
                        Ls = rozdzielone_wsp[5]

                        Ldeg=Ld+Lm/60+Ls/3600
                        L.append(Ldeg)
                        
                        wys = rozdzielone_wsp[6]
                        H.append(wys)
                        
                for (b,l,h) in zip(B,L,H):
                    #utworzy obiekt o tych współrzędnych
                    geodezyjne2000_wgs84=Transformacje(model = "wgs84")
                    X2000, Y2000 = geodezyjne2000_wgs84.BLHto2000(float(b), float(l), float(h))
                    print(f'Wynik transformacji BLH do układu 2000: X2000{X2000:.3f}, Y2000{Y2000:.3f}')
                    wynik_transformacji.write(X2000, Y2000)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
        elif args.t =='BLHto1992':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                B=[]
                L=[]
                H=[]
                for line in lines:
                    if "geodezyjne" in line:
                        continue
                    else:
                        rozdzielone_wsp = line.split(',')
                        Bd = rozdzielone_wsp[0]
                        Bm = rozdzielone_wsp[1]
                        Bs = rozdzielone_wsp[2]
                        
                        Bdeg=Bd+Bm/60+Bs/3600
                        B.append(Bdeg)
                        
                        Ld = rozdzielone_wsp[3]
                        Lm = rozdzielone_wsp[4]
                        Ls = rozdzielone_wsp[5]

                        Ldeg=Ld+Lm/60+Ls/3600
                        L.append(Ldeg)
                        
                        wys = rozdzielone_wsp[6]
                        H.append(wys)
                        
                for (b,l,h) in zip(B,L,H):
                    #utworzy obiekt o tych współrzędnych
                    geodezyjne1992_wgs84=Transformacje(model = "wgs84")
                    X1992, Y1992 = geodezyjne1992_wgs84.BLHto1992(float(b), float(l), float(h))
                    print(f'Wynik transformacji BLH do układu 1992: X1992{X1992:.3f}, Y1992{Y1992:.3f}')
                    wynik_transformacji.write(X1992, Y1992)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
                
    elif args.m == 'krassowski':
        if args.t == 'Hirvonen':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                X=[]
                Y=[]
                Z=[]
                for line in lines:
                    if "geocentryczne" in line:
                        continue
                    else:
                        rozdzielone_wsp=line.split(',')
                        X.append(rozdzielone_wsp[0])
                        Y.append(rozdzielone_wsp[1])
                        Z.append(rozdzielone_wsp[2])
        
                for (x,y,z) in zip(X,Y,Z):
                    #utworzy obiekt o tych współrzędnych
                    geocentryczne_krassowski=Transformacje(model = "krassowski")
                    phi,lam,h = geocentryczne_krassowski.Hirvonen(float(x), float(y), float(z))
                    wynik_transformacji.write(phi,lam,h)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
        elif args.t =='BLHto2000':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                B=[]
                L=[]
                H=[]
                for line in lines:
                    if "geodezyjne" in line:
                        continue
                    else:
                        rozdzielone_wsp = line.split(',')
                        Bd = rozdzielone_wsp[0]
                        Bm = rozdzielone_wsp[1]
                        Bs = rozdzielone_wsp[2]
                        
                        Bdeg=Bd+Bm/60+Bs/3600
                        B.append(Bdeg)
                        
                        Ld = rozdzielone_wsp[3]
                        Lm = rozdzielone_wsp[4]
                        Ls = rozdzielone_wsp[5]

                        Ldeg=Ld+Lm/60+Ls/3600
                        L.append(Ldeg)
                        
                        wys = rozdzielone_wsp[6]
                        H.append(wys)
                        
                for (b,l,h) in zip(B,L,H):
                    #utworzy obiekt o tych współrzędnych
                    geodezyjne2000_krassowski=Transformacje(model = "krassowski")
                    X2000, Y2000 = geodezyjne2000_krassowski.BLHto2000(float(b), float(l), float(h))
                    print(f'Wynik transformacji BLH do układu 2000: X2000{X2000:.3f}, Y2000{Y2000:.3f}')
                    wynik_transformacji.write(X2000, Y2000)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
        elif args.t =='BLHto1992':
            with open(args.from_file, 'r'):
                lines = args.from_file.readlines()
                wynik_transformacji = open(args.to_file, 'w')
                
                B=[]
                L=[]
                H=[]
                for line in lines:
                    if "geodezyjne" in line:
                        continue
                    else:
                        rozdzielone_wsp = line.split(',')
                        Bd = rozdzielone_wsp[0]
                        Bm = rozdzielone_wsp[1]
                        Bs = rozdzielone_wsp[2]
                        
                        Bdeg=Bd+Bm/60+Bs/3600
                        B.append(Bdeg)
                        
                        Ld = rozdzielone_wsp[3]
                        Lm = rozdzielone_wsp[4]
                        Ls = rozdzielone_wsp[5]

                        Ldeg=Ld+Lm/60+Ls/3600
                        L.append(Ldeg)
                        
                        wys = rozdzielone_wsp[6]
                        H.append(wys)
                        
                for (b,l,h) in zip(B,L,H):
                    #utworzy obiekt o tych współrzędnych
                    geodezyjne1992_krassowski=Transformacje(model = "krasowski")
                    X1992, Y1992 = geodezyjne1992_krassowski.BLHto1992(float(b), float(l), float(h))
                    print(f'Wynik transformacji BLH do układu 1992: X1992{X1992:.3f}, Y1992{Y1992:.3f}')
                    wynik_transformacji.write(X1992, Y1992)
                    
                args.from_file.close()
                wynik_transformacji.close()
                
                
            
                