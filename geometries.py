import numpy as np
from utilities import pole, P2R, plot_bb
import svgpathtools

'''
Object hierarchy:

svg input ->
-List_of_Domains, List
    -Domain, Medial axis transform will be calculated for every domain
        -Exterior, List
        -List_of_Interiors, List
            -Loops, List
                -Curves (Line, Arc, Bezier)
'''


'''
Domains
'''
class Domain():
    '''
    A Domain is a single closed area. 
    Its exterior is a closed loop, which can consist of several curve-segments.
    It can also consist of a series of internal holes
    '''
    def __init__(self, exterior, interior=[]):
        self.exterior = exterior
        self.interior = interior

    def __repr__(self):
        return ("Domain")
    
'''
Curves
'''
class Bezier():
    def __init__(self, points):
        assert np.array(points).shape==(4,)
        self.points = np.array(points)
        self.start_point = self.points[0]
        self.end_point = self.points[-1]
        #self.coeffs=[a,b,c,d] with f(t) = a * t^3 + lb * t^2 + c * t^1 + d
        self.coeffs = np.einsum("i,ij",  self.points, np.array([[-1, 3, -3, 1], [3, -6, 3, 0], [-3, 3, 0, 0], [1, 0, 0, 0]]))
    
    def __repr__(self):
        return ("Bezier")

    def find_closest_point(self, point):        
        #candidates for values of t are 0, 1 and all extrema of the distance-function
        candidates=[0,1]
        
        da = 6*np.vdot(self.coeffs[0], self.coeffs[0])
        db = 10*np.vdot(self.coeffs[0], self.coeffs[1])
        dc = 4*(np.vdot(self.coeffs[1], self.coeffs[1]) + 2*np.vdot(self.coeffs[0], self.coeffs[2]))
        dd = 6*(np.vdot(self.coeffs[0], self.coeffs[3]-point) + np.vdot(self.coeffs[1], self.coeffs[2]))
        de = 2*np.vdot(self.coeffs[2], self.coeffs[2]) + 4*np.vdot(self.coeffs[1], self.coeffs[3]-point)
        df = 2*np.vdot(self.coeffs[2], self.coeffs[3]-point)
        
        dcoeffs = np.stack([da, db, dc, dd, de, df]).real
        r = np.roots(dcoeffs)
        extrema = r.real[abs(r.imag)<1e-6]
        for tn in extrema[(extrema>0) & (extrema<1)]:
            candidates.append(tn)
        candidates = np.array(candidates)
        return self.at(candidates[(np.abs(self.at(candidates)-point)).argmin()])
        
    def at(self, t):
        #returns the bezier curve at points 0<t<1
        if isinstance(t, np.ndarray):
            return np.einsum("i,ik->k", self.coeffs, np.stack((t**3,t**2,t**1,t**0)))
        else:
            return np.einsum("i,i", self.coeffs, np.stack((t**3,t**2,t**1,t**0)))
    
    def loose_bounding_box(self):
        # the loose bounding box is definded as the smalest box containing the 4 bezier points        
        self.lbb = np.array([self.points.real.min()+1j*self.points.imag.min(), self.points.real.max()+1j*self.points.imag.max()])
        return self.lbb
    
    def bounding_box(self):
        x_candidates = [self.points[0].real, self.points[-1].real]      
        y_candidates = [self.points[0].imag, self.points[-1].imag]
        
        a = -3*self.points[0]+9*self.points[1]-9*self.points[2]+3*self.points[3]
        b =  6*self.points[0]-12*self.points[1]+6*self.points[2]
        c = -3*self.points[0]+3*self.points[1]

        p_r_roots = np.polynomial.Polynomial([c.real,b.real,a.real]).roots()
        p_i_roots = np.polynomial.Polynomial([c.imag,b.imag,a.imag]).roots()

        for tn in p_r_roots[(p_r_roots>0) & (p_r_roots<1)]:
            x_candidates.append(self.at(tn).real)
        for tn in p_i_roots[(p_i_roots>0) & (p_i_roots<1)]:
            y_candidates.append(self.at(tn).imag)
        
        self.tbb = np.array([min(x_candidates)+1j*min(y_candidates), max(x_candidates)+1j*max(y_candidates)])
        self.tbb_area = np.abs((self.tbb.real[0]-self.tbb.real[1])*(self.tbb.imag[0]-self.tbb.imag[1]))
        return self.tbb
    
    def highest_point(self):
        hp = self.points[0]
        if self.points[-1].imag > hp.imag: hp = self.points[-1]        
    
        a = -3*self.points[0]+9*self.points[1]-9*self.points[2]+3*self.points[3]
        b =  6*self.points[0]-12*self.points[1]+6*self.points[2]
        c = -3*self.points[0]+3*self.points[1]
        p_i_roots = np.polynomial.Polynomial([c.imag,b.imag,a.imag]).roots()
        for tn in p_i_roots[(p_i_roots>0) & (p_i_roots<1)]:
            if self.at(tn).imag > hp.imag: hp = self.at(tn)
        return hp
    
    
class Line():
    def __init__(self, points):
        assert np.array(points).shape==(2,)
        self.points = np.array(points)
        self.start_point = self.points[0]
        self.end_point = self.points[-1]
        #self.coeffs=[a,b] with f(t) = a * t + b
        self.coeffs = np.array([points[1]-points[0], points[0]])

    def __repr__(self):
        return ("Line")
    
    def closest_point(self, point):
        pp = point - self.coeffs[1] 
        t = (self.coeffs[0].real * pp.real + self.coeffs[0].imag * pp.imag)/(self.coeffs[0].real**2 + self.coeffs[0].imag**2)
        if   t<=0: return self.coeffs[1]
        elif t>=1: return self.coeffs[0]+self.coeffs[1]
        else:      return self.at(t)
    
    def at(self, t):
        #returns the bezier curve at points 0<t<1
        return self.coeffs[0] * t + self.coeffs[1]
    
    def bounding_box(self):
        self.tbb = np.array([self.points.real.min()+1j*self.points.imag.min(), self.points.real.max()+1j*self.points.imag.max()])
        self.tbb_area = np.abs((self.tbb.real[0]-self.tbb.real[1])*(self.tbb.imag[0]-self.tbb.imag[1]))
        return self.tbb    
    
    def highest_point(self):
        hp = self.points[0]
        if self.points[-1].imag > hp.imag: hp = self.points[-1]     
        return hp
    
class Arc():
    def __init__(self, start_point, end_point, radius, angle_to_center=0):
        """
        the Arc path leads from start_point to end_point with a constant curvature 1/radius, 
        where radius>0 defines a left bended (counterclockwise) curve and r<0 defines a right bended (clockwise) curve
        For internal calculations only, curves are converted to counterclockwise movement 
        with (s_point, e_point)=(end_point, start_point) if r<0,  else: (s_point, e_point)=(start_point, end_point)
        angles (self.alpha_1, self.alpha_2) refer to (self.s_point, self.e_point)
        """
        #print(end_point, start_point, radius)
        if np.abs(end_point-start_point) > np.abs(2*radius): # this can happen due to rounding
            radius = np.abs(end_point-start_point)/2 * (np.abs(radius)/radius)
        if start_point==end_point:
            self.center = start_point + P2R(radius, angle_to_center)
        else:
            self.center = P2R(np.sqrt(radius**2-np.abs((end_point - start_point)/2)**2), \
                            np.angle(end_point - start_point) + np.pi/2 - np.pi*(radius<0)) + \
                        (end_point + start_point)/2

        self.radius = radius # can be positive or negative
        self.start_point = start_point
        self.end_point = end_point
        
        if self.radius>0:
            self.s_point, self.e_point = start_point, end_point
        else:
            self.s_point, self.e_point = end_point, start_point

        self.alpha_1 = np.angle(self.s_point - self.center)%(2*np.pi) # 0=<alpha_1<2*pi
        self.alpha_2 = np.angle(self.e_point - self.center)%(2*np.pi) # 0 =< alpha_1 < alpha_2 < 4*pi
        self.alpha_2 += 2*np.pi*(self.alpha_1>=self.alpha_2)
            # angles are measured counterclockwise from the positive x-axis

    def __repr__(self):
        return ("Arc")    
    
    def closest_point(self, point):
        alpha_p = np.angle(point - self.center)%(2*np.pi) # 0=<alpha<2*pi        
        alpha_3 = ((self.alpha_1 + self.alpha_2)/2 + np.pi) %(2*np.pi) # 0=<alpha_3<2*pi

        # closest_point lies on the arc
        if (self.alpha_1 <= alpha_p and alpha_p <= self.alpha_2) or (self.alpha_1 <= alpha_p+2*np.pi and alpha_p+2*np.pi <= self.alpha_2):
            if self.radius>0:
                return self.at((alpha_p - self.alpha_1) / (self.alpha_2 - self.alpha_1))
            else:
                return self.at(1- (alpha_p - self.alpha_1) / (self.alpha_2 - self.alpha_1))

        alpha_bool = ((alpha_p - alpha_3)%(2*np.pi))*(self.s_point!=self.e_point)
        
        if alpha_bool > np.pi:       
            return self.e_point
        elif alpha_bool > 0:
            return self.s_point            

         
    def at(self, t):     
        '''
        t=0: self.start_point
        t=1: self.end_point
        '''
        if self.radius>0:    
            return P2R(self.radius, self.alpha_1 * (1-t) + self.alpha_2 * t) + self.center
        else:    
            return P2R(-self.radius, self.alpha_2 * (1-t) + self.alpha_1 * t) + self.center
        
    def bounding_box(self):
        x_candidates = [self.start_point.real, self.end_point.real]      
        y_candidates = [self.start_point.imag, self.end_point.imag]
        a = (self.start_point-self.center).imag>=0
        b = (self.start_point-self.center).real>=0
        c = (self.end_point-self.center).imag>=0
        d = (self.end_point-self.center).real>=0
        if pole(a, b, c, d): 
            y_candidates.append((self.center + 1j*self.radius).imag)
        if pole(not b, a, not d, c): 
            x_candidates.append((self.center - self.radius).real)
        if pole(not a, not b, not c, not d): 
            y_candidates.append((self.center - 1j*self.radius).imag)
        if pole(b, not a, d, not c): 
            x_candidates.append((self.center + self.radius).real)
        if self.start_point == self.end_point:
            x_candidates.append((self.center - self.radius).real)
            x_candidates.append((self.center + self.radius).real)
            y_candidates.append((self.center - 1j*self.radius).imag)
            y_candidates.append((self.center + 1j*self.radius).imag)
            
        self.tbb = np.array([min(x_candidates)+1j*min(y_candidates), max(x_candidates)+1j*max(y_candidates)])
        self.tbb_area = np.abs((self.tbb.real[0]-self.tbb.real[1])*(self.tbb.imag[0]-self.tbb.imag[1]))
        return self.tbb
    
    def highest_point(self):
        if self.start_point == self.end_point: return self.center + 1j*self.radius
        hp = self.start_point
        if self.end_point.imag > hp.imag: hp = self.end_point

        a = (self.start_point-self.center).imag>=0
        b = (self.start_point-self.center).real>=0
        c = (self.end_point-self.center).imag>=0
        d = (self.end_point-self.center).real>=0
        if pole(a, b, c, d): 
            if (self.center + 1j*self.radius).imag > hp.imag: hp = self.center + 1j*self.radius
        return hp
    
'''
functions
'''

def bb1_cut_bb2(bb1, bb2):
    # returns True if the two bounding boxes intersect
    return not (bb1[1].real<bb2[0].real or bb2[1].real<bb1[0].real or bb1.imag[1]<bb2.imag[0] or bb2.imag[1]<bb1.imag[0])

def split_half(curve):
    if isinstance(curve, Bezier):
        p0, p1, p2, p3 = curve.points
        p4 = 0.5*p0 + 0.5*p1
        p5 = 0.5*p1 + 0.5*p2
        p6 = 0.5*p2 + 0.5*p3
        p7 = 0.5*p4 + 0.5*p5
        p8 = 0.5*p5 + 0.5*p6
        p9 = 0.5*p7 + 0.5*p8
        return Bezier([p0, p4, p7, p9]), Bezier([p9, p8, p6, p3])
    
    elif isinstance(curve, Line):
        p0, p2 = curve.points
        p1 = (p0+p2)/2
        return Line([p0, p1]), Line([p1, p2])    
    
    elif isinstance(curve, Arc):        
        #print("curve.start_point, curve.at(0.5), curve.radius: ", curve.start_point, curve.at(0.5), curve.radius)
        return Arc(curve.start_point, curve.at(0.5), curve.radius), Arc(curve.at(0.5), curve.end_point, curve.radius)
    
def find_intersection(obj1, obj2):    
    threshold=1e-10
    if not bb1_cut_bb2(obj1.bounding_box(), obj2.bounding_box()): 
        return []
    circumference = (obj1.tbb[1]-obj1.tbb[0]).real+(obj1.tbb[1]-obj1.tbb[0]).imag + \
        (obj2.tbb[1]-obj2.tbb[0]).real+(obj2.tbb[1]-obj2.tbb[0]).imag
    if circumference < threshold:        
        plot_bb(obj1.bounding_box())        
        plot_bb(obj2.bounding_box())
        return [(obj1.tbb[0]+obj1.tbb[1])/2]
    obj1a, obj1b = split_half(obj1)
    obj2a, obj2b = split_half(obj2)
    a, b, c, d = find_intersection(obj1a, obj2a), find_intersection(obj1a, obj2b), find_intersection(obj1b, obj2a), find_intersection(obj1b, obj2b)
    intersections_found = []
    n=0
    if a: 
        intersections_found += a  
        n+=1
    if b: 
        intersections_found += b 
        n+=1 
    if c: 
        intersections_found += c  
        n+=1
    if d: 
        intersections_found += d  
        n+=1
    if n>=2:        
        del_list=[]
        for p1i in range(len(intersections_found)):
            for p2i in range(p1i+1, len(intersections_found)):
                if np.abs(intersections_found[p1i]-intersections_found[p2i]) < threshold*2:
                    del_list.append(p2i)
        del_list = sorted(list(set(del_list)))
        for pi in del_list[::-1]:
            del intersections_found[pi]
    return intersections_found

def is_hole(closed_loop):
    soa = sum_of_angles(closed_loop)
    return soa==np.abs(soa)

def sum_of_angles(closed_loop):
    #print("len(closed_loop): ", len(closed_loop))
    angle_sum=0
    angle = np.angle(closed_loop[0].end_point-closed_loop[0].start_point)-np.angle(closed_loop[-1].end_point-closed_loop[-1].start_point)
    angle_sum+=(angle+np.pi)%(2*np.pi)-np.pi   
    for n in range(1, len(closed_loop)):
        angle = np.angle(closed_loop[n].end_point-closed_loop[n].start_point)-np.angle(closed_loop[n-1].end_point-closed_loop[n-1].start_point)
        angle_sum+=(angle+np.pi)%(2*np.pi)-np.pi
    return angle_sum

def create_list_of_domains(svgPaths):
    list_of_domains = []
    for dn in range(len(svgPaths)):
        exterior = []
        interior = []
        loop = [] 
        start_of_loop = svgPaths[dn]._segments[0].start
        cursor = start_of_loop
        for n, s in enumerate(svgPaths[dn]._segments):
            #print(type(s))
            if isinstance(s, svgpathtools.path.Line):
                loop.append(Line([s.start, s.end]))
                cursor = s.end            
            elif isinstance(s, svgpathtools.path.CubicBezier):
                loop.append(Bezier([s.start, s.control1, s.control2, s.end]))
                cursor = s.end
            elif isinstance(s, svgpathtools.path.QuadraticBezier):
                loop.append(Bezier([s.start, s.start+2/3*(s.control-s.start), s.end+2/3*(s.control-s.end), s.end]))
                cursor = s.end
            elif isinstance(s, svgpathtools.path.Arc):
                print("Input of Arcs not yet implemented")

            if cursor == start_of_loop:
                # print(loop, is_hole(loop))
                if is_hole(loop): interior.append(loop)
                else: exterior=loop
                if n!=len(svgPaths[dn]._segments)-1:
                    loop=[]
                    start_of_loop = svgPaths[dn]._segments[n+1].start
        if len(exterior)>0:
            list_of_domains.append(Domain(exterior, interior))
    return list_of_domains

def print_domain_structure(list_of_domains):
    for D in list_of_domains:
        print(D)
        print("    Exterior: ")
        print("    ", D.exterior)
        print("    Interior: ")
        for I in D.interior:        
            print("    ", I)
    return None

