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
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} niezaimplementowany model")
        self.flat = (self.a - self.b) / self.a
        self.ecc = np.sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
        
        
    def Hirvonen(self, X, Y, Z, output = 'dec_degree'):
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
        phi_poprzednie= np.arctan(Z / (p * (1 - self.ecc2)))    # pierwsze przybliilizenie
        phi = 0
        while abs(phi_poprzednie - phi) > 0.000001/206265:    
            phi_poprzednie = phi
            N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi_poprzednie)**2)
            h = p / np.cos(phi_poprzednie) - N
            phi = np.arctan((Z/p) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lam = np.arctan(Y/X)
        N = self.a / np.sqrt(1 - self.ecc2 * (np.sin(phi))**2);
        h = p / np.cos(phi) - N       
        if output == "dec_degree":
            def deg2dms(dd):
                deg=np.trunc(dd)
                mnt=np.trunc((dd-deg)*60) 
                sec=(((dd-deg)*60)-mnt)*60
                print(f'{deg:.0f} {abs(mnt):.0f} {abs(sec):.5f}')
            return np.rad2deg(phi), np.rad2deg(lam), h 
    
        elif output == "dms":
            phi = self.deg2dms(np.rad2deg(phi))
            lam = self.deg2dms(np.rad2deg(lam))
            return f"{phi[0]:02d}:{phi[1]:02d}:{phi[2]:.2f}", f"{lam[0]:02d}:{lam[1]:02d}:{lam[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - format niezdefiniowany")
            
    
            
            
    def BLHto2000(self, phi, lam, h):
            '''
            
            Transformacja współrzędnych geodezyjnych φ i λ n
            a brytyjski układ współrzędnych geodezyjnych UK2000 
            to przekształcenie z jednego układu odniesienia na drugi. 
            Wartości φ i λ wyrażają położenie punktu na powierzchni elipsoidy referencyjnej, a układ UK2000 jest planimetrycznym układem współrzędnych geodezyjnych, który stosuje się do map, a jego wartości współrzędnych wyrażone są w metrach.

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
        
    def XYZ2neu(self, xa, ya, za, xb, yb, zb, phi, lam, h):
        '''
        XYZ -> NEUp -
        Transformacja z układu współrzędnych XYZ (układu kartezjańskiego) 
        na układ współrzędnych NEU (układ lokalny związany z punktem odniesienia)
        służy do dokładnego okreslenia położenia punktu w lokalnym układzie współrzędnych.



        '''
        # współrzędne początku  φ λ -> XYZ 
        N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
        X0 = (N + h) * np.cos(phi) * np.cos(lam)
        Y0 = (N + h) * np.cos(phi) * np.sin(lam)
        Z0 = ((1 - self.ecc2) * N + h) * np.sin(phi)
    
        # obliczenie współrzędnych XYZ
        dxyz = np.array([Xb,Yb,Zb]) - np.array([Xa,Ya,Za])
        X, Y, Z = dxyz - np.array([X0, Y0, Z0])
    
        # Macierz obrotu
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        sin_lam = np.sin(lam)
        cos_lam = np.cos(lam)
    
        R = np.array([[-sin_lam,           cos_lam,           0],
                      [-sin_phi*cos_lam, -sin_phi*sin_lam,  cos_phi],
                      [ cos_phi*cos_lam,  cos_phi*sin_lam,  sin_phi]])
    
        # Obliczenie współrzędnych NEU
        NEU = np.dot(R, np.array([X, Y, Z]))
    
     
        N = NEU[0]
        E = NEU[1]
        U = -NEU[2]
    
    
        phi_stopnie = np.degrees(phi)
        lam_stopnie = np.degrees(lam)
    
        return (N, E, U, phi_stopnie, lam_stopnie)

if __name__=='__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='Transformacje wspolrzednych')
    parser.add_argument('x', type=float, help='Współrzędna x')
    parser.add_argument('y', type=float, help='Współrzędna y')
    parser.add_argument('z', type=float, help='Współrzędna z')
    args = parser.parse_args()
    
    phi,lam,h = 
    
    
    with open('plikprzykladowedane.py', 'r') as f:
        lines = f.readlines()

        X = []
        Y = []
        Z = []
        
        PHI=[]
        LAM=[]
        H=[]
        
        N=[]
        E=[]
        U=[]
        
        X2000=[]
        Y2000=[]

        
        X1992=[]
        Y1992=[]
        
        
        for line in lines:
            if "print" in line:
                continue
            else:
                rozdzielone_wsp=line.split(',')
                X.append(rozdzielone_wsp[0])
                Y.append(rozdzielone_wsp[1])
                Z.append(rozdzielone_wsp[2])
        
        for (x,y,z) in zip(X,Y,Z):
        #utworzy obiekt o tych współrzędnych
            geocentryczne=Transformacje(model = "grs80")
            phi,lam,h = geocentryczne.Hirvonen(float(x), float(y), float(z))
            
            print('Wynikiem transformacji Hirvonena dla podanych X,Y,Z są współrzędne: ')
            print(f'B={phi} stopnie, L={lam}stopnie, H={h}m ')
            
            geodezyjne=Transformacje(model = "wgs84")
            X2000,Y2000=geodezyjne.BLHto2000(phi, lam, h)
            print('Wynikiem transformacji BLH do układu 2000 dla podanych phi, lam, h są współrzędne: ')
            print(f'X2000={X2000:.3f}m , Y2000={Y2000:.3f}m ')
            
            geodezyjne = Transformacje(model = "wgs84")
            X1992,Y1992 = geodezyjne.BLHto1992(phi, lam, h)
            print('Wynikiem transformacji BLH do układu 1992 dla podanych phi, lam, h są współrzędne: ')
            print(f'X1992={X1992:.3f}m , Y1992={Y2000:.3f}m ')
            print("=================")
            print("=================")