from wtforms import StringField,SubmitField,PasswordField
from wtforms.validators import  DataRequired,Length
from flask_wtf import FlaskForm

#登录表单
class Login_Form(FlaskForm):
    # name=StringField('name',validators=[ DataRequired(u"用户名不能为空")])
    # pwd=PasswordField('pwd',validators=[DataRequired(u"密码不能为空"), Length(6,16,message='密码长度应在6-16位')])
    submit=SubmitField(label='Login in')

#注册表单
class Register_Form(FlaskForm):
    name=StringField('name',validators=[DataRequired(u"用户名不能为空")])
    pwd=PasswordField('pwd',validators=[DataRequired(u"密码不能为空"), Length(6,16,message='密码长度应在6-16位')])
    submit=SubmitField('register')


# #导入表单基类
# from flask_wtf import FlaskForm
# #导入需要用到的表格类型
# from wtforms import StringField,PasswordField,SubmitField
# #导入需要用到的验证方式
# from wtforms.validators import DataRequired,EqualTo,Length
#
# #用户注册表单
# class RegisterForm(FlaskForm):
#     #用户名表框
#     username = StringField(
#         label= '用户名',    # 标签
#         validators=[        # 验证方式
#             DataRequired(u"用户名不能为空")  # 不能为空
#         ]
#     )
#     # 密码表框
#     password = PasswordField(
#         label='密码',
#         validators=[
#             DataRequired(u"密码不能为空"),
#             Length(6,16,message='密码长度应在6-16位') # 密码长度必须在6-16，否则输出提示信息
#         ]
#     )
#     # 注册按钮
#     submit = SubmitField(
#         label='注册'
#     )
#
# #登录表单
# class LoginForm(FlaskForm):
#     username=StringField(
#         label='用户名',
#         validators=[DataRequired(u"用户名不能为空")]
#     )
#
#     password=PasswordField(
#         label='密码',
#         validators=[DataRequired(u"密码不能为空"),
#             Length(6,16,message='密码长度应在6-16位')]
#     )
#
#     submit=SubmitField(
#         label='登录'
#     )
#
#
#
#
#
#
#
#
#
#
#
#
