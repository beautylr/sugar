#导入表单基类
from flask_wtf import FlaskForm
#导入需要用到的表格类型
from wtforms import SubmitField
#导入需要用到的验证方式
from wtforms.validators import DataRequired,EqualTo,Length


class AppForm(FlaskForm):
    submit=SubmitField(
        label='提交评价',
    )

class Trajectory(FlaskForm):
    submit=SubmitField(
        label='行驶轨迹',
    )

class Error(FlaskForm):
    submit=SubmitField(
        label='误差分析',
    )







