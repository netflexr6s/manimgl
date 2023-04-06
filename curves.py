from manimlib import *
from math import cos, sin

# Parametrisierungen
circ = lambda t: np.array([cos(t),sin(t),0])
doublecirc = lambda t: np.array([cos(2*t),sin(2*t),0])
lemniscate = lambda t: np.array([3*cos(t) / (1+sin(t)**2), 3*sin(t)*cos(t) / (1+sin(t)**2),0])


class TracedCurve(Scene):
    def construct(self):
        ax = NumberPlane()
        dot = GlowDot()
        
        self.func = lemniscate
        self.range = (0,2*PI)
        
        curve = ParametricCurve(
            self.func,
            self.range
        )

        parameter_space = NumberLine(
            (0,np.ceil(self.range[1])),
            include_numbers = True
        ).set_color(WHITE).scale(0.75)
        parameter_space.to_corner(DR)

        parameter_dot = Dot().set_color(YELLOW).move_to(parameter_space.n2p(0))
        parameter = DecimalNumber(
            0, font_size = 38, 
            num_decimal_places = 2
        )
        parameter.add_updater(
            lambda x: x.set_value(
                parameter_space.p2n(
                    parameter_dot.get_center()
                )
            )
        )
        t_label = Tex("t=", font_size = 38)
        par_label = VGroup(t_label, parameter).set_color(YELLOW)
        par_label.arrange(RIGHT, buff=SMALL_BUFF)
        always(par_label.next_to, parameter_dot, UP)

        parameter_path = Line(
            parameter_space.n2p(0),
            parameter_space.n2p(self.range[1])
        )
        par = VGroup(
            parameter_space, parameter_dot, par_label
        )

        trace = TracedPath(
            lambda: dot.get_center(), 
            stroke_width = 3
        ).set_color(RED)

        dot.move_to(self.func(0))
        self.add(ax, trace, dot)
        self.add(par)
        self.interact()

        self.play(
            MoveAlongPath(
                dot, curve
            ),
            MoveAlongPath(
                parameter_dot, parameter_path
            ),
            run_time = 5
        )


    # Interaktion
        dot2 = GlowDot()
        dot2.add_updater(
            lambda x: x.move_to(
                self.func(
                    parameter_space.p2n(
                        parameter_dot.get_center()
                    )
                )
            )
        )
        self.remove(dot)
        self.add(dot2)

        self.search_set = VGroup(parameter_dot)
        self.line = parameter_space

    def on_mouse_press(self, point, button, mods):
        super().on_mouse_press(point, button, mods)
        mob = self.point_to_mobject(point, search_set=self.search_set)
        if mob is None:
            return
        self.mouse_drag_point.move_to(point) 
        mob.add_updater(
            lambda x: x.set_x(
                self.mouse_drag_point.get_center()[0]
            ) 
        )
        self.unlock_mobject_data()
        self.lock_static_mobject_data()

    def on_mouse_release(self, point, button, mods):
        super().on_mouse_release(point, button, mods)
        self.search_set.clear_updaters() 



# Parametrisierungen
helix = lambda t: np.array([cos(t),sin(t),np.sqrt(t)])
spherical1 = lambda t: np.array([sin(2*t)*cos(3*t),sin(2*t)*sin(3*t),cos(2*t)])
spherical2 = lambda t: np.array([sin(3*t)*cos(5*t),sin(3*t)*sin(5*t),cos(3*t)])


class TracedCurve3D(Scene):
    CONFIG = {
        "camera_class": ThreeDCamera
    }
    def construct(self):
        plane = NumberPlane(
            (-3,3),
            (-3,3),
            width = 6,
            height = 6
        ).set_opacity(0.5)
        ax = ThreeDAxes(
            (-3,3),
            (-3,3),
            (-3,3),
            width = 6,
            height = 6
        )
        dot = GlowDot()

        self.func = spherical2
        self.range = (0,2*PI)

        trace = TracedPath(
            lambda: dot.get_center(), 
            stroke_width = 3
        ).set_color(RED)

        curve = ParametricCurve(
            self.func,
            self.range
        )

        parameter_space = NumberLine(
            (0,np.ceil(self.range[1])),
            include_numbers = True
        ).set_color(WHITE).scale(0.75)
        parameter_space.to_corner(DR)

        parameter_dot = Dot().set_color(YELLOW).move_to(parameter_space.n2p(0))
        parameter = DecimalNumber(
            0, font_size = 38, 
            num_decimal_places = 2
        )
        parameter.add_updater(
            lambda x: x.set_value(
                parameter_space.p2n(
                    parameter_dot.get_center()
                )
            )
        )
        t_label = Tex("t=", font_size = 38)
        par_label = VGroup(t_label, parameter).set_color(YELLOW)
        par_label.arrange(RIGHT, buff=SMALL_BUFF)
        always(par_label.next_to, parameter_dot, UP)

        parameter_path = Line(
            parameter_space.n2p(0),
            parameter_space.n2p(self.range[1])
        )
        par = VGroup(
            parameter_space, parameter_dot, par_label
        ).fix_in_frame()


        dot.move_to(self.func(0))
        self.add(plane, ax, dot, trace)
        self.add(par)
        self.interact()

        self.play(
            MoveAlongPath(
                dot, curve
            ),
            MoveAlongPath(
                parameter_dot, parameter_path
            ),
            run_time = 5
        )

    # Interaktion
        dot2 = GlowDot()
        dot2.add_updater(
            lambda x: x.move_to(
                self.func(
                    parameter_space.p2n(
                        parameter_dot.get_center()
                    )
                )
            )
        )
        self.remove(dot)
        self.add(dot2)

        self.search_set = VGroup(parameter_dot)
        self.line = parameter_space
        
    def on_mouse_press(self, point, button, mods):
        super().on_mouse_press(point, button, mods)
        mob = self.point_to_mobject(point, search_set=self.search_set)
        if mob is None:
            return
        self.mouse_drag_point.move_to(point) 
        mob.add_updater(
            lambda x: x.set_x(
                self.mouse_drag_point.get_center()[0]
            ) 
        )
        self.unlock_mobject_data()
        self.lock_static_mobject_data()

    def on_mouse_release(self, point, button, mods):
        super().on_mouse_release(point, button, mods)
        self.search_set.clear_updaters()

        

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