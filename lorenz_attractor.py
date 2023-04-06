from manimlib import *
from scipy.integrate import odeint

class LorenzAttractor(Scene):
    CONFIG = {
        "camera_class": ThreeDCamera
    }
    def construct(self):
        frame = self.camera.frame
        frame.set_euler_angles(phi=80*DEGREES,theta=35*DEGREES)  
        frame.scale(0.9)

        def lorenz(t, state, s=10, r=28, b=2.667):
            x, y, z = state
            x_dot = s*(y - x)
            y_dot = r*x - y - x*z
            z_dot = x*y - b*z
            return [x_dot, y_dot, z_dot]
        
        def lorenz_curve(init): 
            grid = np.arange(0.0,45.0,0.01)
            res = odeint(lorenz, init, grid, tfirst=True)
            p = VMobject()
            p.set_points_as_corners([*[[a,b,c] for a,b,c in zip(res[:,0],res[:,1],res[:,2])]])
            p.make_smooth().set_stroke(None,1)
            p.scale(0.1).move_to(ORIGIN)
            return p
        
        colors = [RED, ORANGE, YELLOW] 

        dot1 = GlowDot()
        dot2 = GlowDot()

        trace_tail1 = TracingTail(
            lambda: dot1.get_center(), 
            stroke_width = 3,
            time_traced = 2
        ).set_color_by_gradient(RED)

        trace_tail2 = TracingTail(
            lambda: dot2.get_center(), 
            stroke_width = 3,
            time_traced = 2
        ).set_color_by_gradient(BLUE)

        trace1 = TracedPath(
            lambda: dot1.get_center(), 
            stroke_width = 2,
        ).set_color_by_gradient(RED)

        trace2 = TracedPath(
            lambda: dot2.get_center(), 
            stroke_width = 2,
        ).set_color_by_gradient(BLUE)

        
        #self.add(trace_tail1, trace_tail2)
        self.add(trace1, trace2)
        self.add(dot1, dot2)
        self.play(
            MoveAlongPath(
                dot1, lorenz_curve([0.0,1.0,1.05]),
                rate_func = linear
            ),
            MoveAlongPath(
                dot2, lorenz_curve([0.1,1.1,1.08]),
                rate_func = linear
            ),
            run_time = 40
        )
        self.wait(2)
