# ManimGL v1.6.1
from manimlib import *
from scipy.integrate import odeint
from math import sin, cos

    
class DoublePendulum(VGroup):
    CONFIG = {
        "L1": 2,
        "L2": 1.5,
        "m1": 4,
        "m2": 5,
        "g": 9.8,
        "init_t1": 60*DEGREES,
        "init_t2": 30*DEGREES,
        "init_o1": 2,
        "init_o2": 2,
        "trace_color": RED
    }
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        L1 = self.L1
        L2 = self.L2

        self.dot1 = dot1 = Dot(L1*np.array([sin(self.init_t1),-cos(self.init_t1),0])).set_color(YELLOW)
        self.dot2 = dot2 = Dot(
            L1*np.array([sin(self.init_t1),-cos(self.init_t1),0]) + L2*np.array([sin(self.init_t2),-cos(self.init_t2),0])
        ).set_color(RED)

        line1 = always_redraw(
            lambda: Line(np.array([0,0,0]), dot1.get_center()).set_color(WHITE)
        )
        line2 = always_redraw(
            lambda: Line(dot1.get_center(), dot2.get_center()).set_color(WHITE)
        )

        self.add(line1, line2, dot1, dot2)
        trace = TracingTail(dot2, time_traced=2).set_stroke(self.trace_color, 1) 
        self.add(trace)

    def get_positions(self):
        L1 = self.L1
        L2 = self.L2
        p1 = self.get_angles([self.init_t1, self.init_o1, self.init_t2, self.init_o2])[:,0] 
        p2 = self.get_angles([self.init_t1, self.init_o1, self.init_t2, self.init_o2])[:,2]

        x1 = [L1*sin(t1) for t1 in p1]
        y1 = [-L1*cos(t1) for t1 in p1]

        x2 = [L2*sin(t2)+L1*sin(t1) for t1, t2 in zip(p1, p2)]
        y2 = [-L2*cos(t2)-L1*cos(t1) for t1, t2 in zip(p1, p2)]

        pos1 = [x1[k+1]*RIGHT + y1[k+1]*UP for k in range(len(x1)-1)]
        pos2 = [x2[k+1]*RIGHT + y2[k+1]*UP for k in range(len(x1)-1)]
        return [
            [p1, p2] for p1,p2 in zip(pos1,pos2)
        ]

    
    def get_angles(self, pos):
        L1 = self.L1
        L2 = self.L2
        M1 = self.m1
        M2 = self.m2
        G = self.g
        
        def DGL(t, state):
            dydx = np.zeros_like(state)
            dydx[0] = state[1]
            delta = state[2] - state[0]
            den1 = (M1+M2) * L1 - M2 * L1 * cos(delta) * cos(delta)
            dydx[1] = ((M2 * L1 * state[1] * state[1] * sin(delta) * cos(delta)
                        + M2 * G * sin(state[2]) * cos(delta)
                        + M2 * L2 * state[3] * state[3] * sin(delta)
                        - (M1+M2) * G * sin(state[0]))
                       / den1)
            dydx[2] = state[3]
            den2 = (L2/L1) * den1
            dydx[3] = ((- M2 * L2 * state[3] * state[3] * sin(delta) * cos(delta)
                        + (M1+M2) * G * sin(state[0]) * cos(delta)
                        - (M1+M2) * L1 * state[1] * state[1] * sin(delta)
                        - (M1+M2) * G * sin(state[2]))
                        / den2)
            return dydx

        def path(init): 
            grid = np.arange(0.0, 40, 0.08)
            res = odeint(DGL, init, grid, tfirst=True)
            return res
        
        return path(pos)
    
class DoublePendulumAnimation(Scene):
    def construct(self):
        bg = NumberPlane()
        self.add(bg)
        pend1 = DoublePendulum(
            init_t1 = 95*DEGREES, 
            init_t2 = 30*DEGREES
        )
        pend2 = DoublePendulum(
            init_t1 = 90*DEGREES,
            init_t2 = 35*DEGREES, 
            trace_color = BLUE_E
        )
        self.add(pend1, pend2)
        pos1 = pend1.get_positions()
        pos2 = pend2.get_positions()
        self.wait()
        for p,q in zip(pos1,pos2):
            self.play(
                pend1.dot1.animate.move_to(p[0]),
                pend1.dot2.animate.move_to(p[1]),  
                pend2.dot1.animate.move_to(q[0]),
                pend2.dot2.animate.move_to(q[1]),
                run_time = 0.1,
                rate_func = linear
            )
