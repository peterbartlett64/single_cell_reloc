from manim import *

class ManimCELogo(Scene):
    def construct(self):
        self.camera.background_color = "#ece6e2"
        logo_green = "#87c2a5"
        logo_blue = "#525893"
        logo_red = "#e07a5f"
        logo_black = "#343434"
        ds_m = MathTex(r"\mathbb{M}", fill_color=logo_black).scale(7)
        ds_m.shift(2.25 * LEFT + 1.5 * UP)
        circle = Circle(color=logo_green, fill_opacity=1).shift(LEFT)
        square = Square(color=logo_blue, fill_opacity=1).shift(UP)
        triangle = Triangle(color=logo_red, fill_opacity=1).shift(RIGHT)
        logo = VGroup(triangle, square, circle, ds_m)  # order matters
        logo.move_to(ORIGIN)
        self.add(logo)


# from manim import *
# import numpy as np

# class PercentileAverage(Scene):
#     def construct(self):
#         data = np.random.normal(loc=0, scale=1, size=1000)  # Generate random data

#         # Calculate the 99th percentile value
#         percentile_99 = np.percentile(data, 99)

#         # Filter values greater than or equal to the 99th percentile
#         filtered_data = data[data >= percentile_99]

#         # Calculate the average of the filtered data
#         average = np.median(filtered_data)

#         # Display the data and the average
#         data_text = Text(f"Data: {data}", font_size=20)
#         percentile_text = Text(f"99th Percentile: {percentile_99:.2f}", font_size=20)
#         average_text = Text(f"Average of 99th Percentile: {average:.2f}", font_size=20)

#         self.play(Write(data_text))
#         self.wait()
#         self.play(Write(percentile_text))
#         self.wait()
#         self.play(Write(average_text))
#         self.wait(3)

# PercentileAverage