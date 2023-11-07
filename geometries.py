import numpy as np
from utilities import pole, P2R, plot_bb

class Bezier():
    def __init__(self, points):
        assert np.array(points).shape==(4,)
        self.points = np.array(points)
        #self.coeffs=[a,b,c,d] with f(t) = a * t^3 + lb * t^2 + c * t^1 + d
        self.coeffs = np.einsum("i,ij",  self.points, np.array([[-1, 3, -3, 1], [3, -6, 3, 0], [-3, 3, 0, 0], [1, 0, 0, 0]]))
    
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
    
    def tight_bounding_box(self):
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
    
class Line():
    def __init__(self, points):
        assert np.array(points).shape==(2,)
        self.points = np.array(points)
        #self.coeffs=[a,b] with f(t) = a * t + b
        self.coeffs = np.array([points[1]-points[0], points[0]])
    
    def at(self, t):
        #returns the bezier curve at points 0<t<1
        return self.coeffs[0] * t + self.coeffs[1]
    
    def tight_bounding_box(self):
        self.tbb = np.array([self.points.real.min()+1j*self.points.imag.min(), self.points.real.max()+1j*self.points.imag.max()])
        self.tbb_area = np.abs((self.tbb.real[0]-self.tbb.real[1])*(self.tbb.imag[0]-self.tbb.imag[1]))
        return self.tbb    
    
class Arc():
    def __init__(self, start_point, end_point, radius):
        """
        input:  center (point coords in x + iy)
                start_point (point coords in x + iy)
                end_point (point coords in x + iy)
        """
        assert np.abs(end_point-start_point) <= 2*radius
        self.center = P2R(np.sqrt(radius**2-np.abs((end_point - start_point)/2)**2), \
                   np.angle(end_point - start_point) + np.pi/2) + (end_point + start_point)/2
        self.radius = radius
        self.start_point = start_point
        self.end_point = end_point
        self.alpha_1 = np.angle(self.start_point - self.center) # 0=<alpha_1<2*pi
        self.alpha_2 = np.angle(self.end_point - self.center) + 2*np.pi*(start_point==end_point) # 0 =< alpha_1 < alpha_2 < 4*pi
        self.alpha_2 += 2*np.pi*(self.alpha_1>self.alpha_2)
        # angles counterclockwise angle measured from the positive x-axis
        
    def closest_point(self, point):
        alpha_p = np.angle(point - self.center) # 0=<alpha<2*pi        
        alpha_3 = ((self.alpha_1 + self.alpha_2)/2 + np.pi) % (2 * np.pi) # 0=<alpha_3<2*pi
        
        if self.alpha_1 <= alpha_p and alpha_p <= self.alpha_2:
            return self.at((alpha_p - self.alpha_1) / (self.alpha_2 - self.alpha_1))
        
        alpha_bool = ((alpha_p - alpha_3)%(2*np.pi))*(self.start_point!=self.end_point)
        
        if alpha_bool > np.pi:       
            return self.end_point
        elif alpha_bool > 0:
            return self.start_point            
        else:
            return self.at((alpha_p - self.alpha_1) / (self.alpha_2 - self.alpha_1))
         
    def at(self, t):        
        return P2R(self.radius, self.alpha_1 * (1-t) + self.alpha_2 * t) + self.center
        
    def tight_bounding_box(self):
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
        return Arc(curve.start_point, curve.at(0.5), curve.radius), Arc(curve.at(0.5), curve.end_point, curve.radius)
    
def find_intersection(obj1, obj2):    
    threshold=1e-10
    if not bb1_cut_bb2(obj1.tight_bounding_box(), obj2.tight_bounding_box()): 
        return []
    circumference = (obj1.tbb[1]-obj1.tbb[0]).real+(obj1.tbb[1]-obj1.tbb[0]).imag + \
        (obj2.tbb[1]-obj2.tbb[0]).real+(obj2.tbb[1]-obj2.tbb[0]).imag
    if circumference < threshold:        
        plot_bb(obj1.tight_bounding_box())        
        plot_bb(obj2.tight_bounding_box())
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