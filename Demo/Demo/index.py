#-*- codeing = utf-8 -*-
#@Time :  14:27
#@AUthor : LR
#@File : index.py
#@Software:PyCharm

from flask import Flask  #从flask中导入Flask，用于创建主app
from flask import render_template,redirect,url_for,session,flash,Response
from ca1 import VideoCapture
from flask import request
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import glob
import plotly as py
import plotly.graph_objs as go
import time
from flask_mysqldb import MySQL
from wtforms import Form, StringField, TextAreaField, PasswordField, validators
from functools import wraps
from app_form import Trajectory
from app_form import Error
from app_form import AppForm

# plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False

#Flask(__name__).run()
app = Flask(__name__)  #创建Flask的主app
app.config["SECRET_KEY"]="lH4WHi5amT0ZqykHvLofllRJu3UN1uzmeUN0z2IiacjDUb5TLU3ZTtUP5VJqgkMY"
app.config['MYSQL_HOST'] = 'localhost'
app.config['MYSQL_USER'] = 'root'
app.config['MYSQL_PASSWORD'] = 'lr345678'
app.config['MYSQL_DB'] = 'testdb'
app.config['MYSQL_CURSORCLASS'] = 'DictCursor'
app.config['MYSQL_PORT'] = 3306
# init MYSQL
mysql = MySQL(app)


# 相机推流
def gen(camera):
    while True:
        _,frame = camera.read()
        yield (b'--frame\r\n'
               b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n\r\n')


# 相机喂流
@app.route('/video_feed')
def video_feed():
    cam = "rtsp://admin:kcy123456@172.20.20.102:554/h264/ch1/main/av_stream"
    return Response(gen(VideoCapture(cam)),
                    mimetype='multipart/x-mixed-replace; boundary=frame')

#系统进入首页
@app.route('/',methods=['GET','POST'])
def trajectory_none():

    return render_template('trajectory_none.html')

#视频监控模块
@app.route('/video_monitor',methods=['GET','POST'])

def video_monitor():

    return render_template('video_monitor.html')

#测试场地监控模块
@app.route('/map_video',methods=['GET','POST'])

def map_video():


    return render_template('map_video_test.html')

@app.route('/scene1',methods=['GET','POST'])
#场景模块
def scene1():
    return render_template('scene1.html')

@app.route('/scene2',methods=['GET','POST'])
#场景模块
def scene2():
    return render_template('scene2.html')

@app.route('/scene3',methods=['GET','POST'])


@app.before_first_request
def before_first_request():
    session['click'] = -1

@app.route('/get_data',methods=['GET','POST'])
#读取OBU数据模块
def obu_data():
    # time.sleep(5)
    # time_c,lon,lat,height,speed,angle,vehtype,light,brake,maxspeed,vehstate,printinfor = testdis()
    time_c,vehtype,vehstate,speed,lon,lat,height,angle,light,brake,maxspeed,printinfor = testdis()
    stop = 0
    # if s == 'z':
    #     stop = 1
    return render_template("obu_data.html",time_c = time_c,lon=lon,lat=lat,height=height,speed=speed,
            angle=angle,vehtype=vehtype,light=light,brake=brake,maxspeed=maxspeed,vehstate=vehstate,printinfor=printinfor,stop=stop)

def testdis():
    filename = "a_425.log"
    session['click']+=1
    print(session['click'])
    # data=linecache.getline(filename,session['click'])
    # return data
    with open(filename, 'r') as csvfile:
        # reader = csv.reader(csvfile)
        # reader=csvfile.readline()
        # print(reader)
        for i, rows in enumerate(csvfile):
            if i==session['click']:
                if rows.split("\t")[0].strip()=="DATA_FLAG".strip():
                    str="OBU Data Parsed Success!"
                    rows=rows.split("\t")
                    return rows[1],rows[2],rows[3],rows[4],rows[5],rows[6],rows[7],rows[8],rows[9],rows[10],rows[11],str
                else:
                    return 0,0,0,0,0,0,0,0,0,0,0,rows


@app.route('/map',methods=['GET','POST'])
#地图显示
def map():

    return render_template("map.html",title='Home')

#车辆行驶参数模块-行驶轨迹
@app.route('/trajectory',methods=['GET','POST'])
def trajectory():

    return render_template("trajectory.html",title='Home') #plotly画图
    # return render_template("trajectory_echarts.html",human_trace_x=human_trace_x,human_trace_y=human_trace_y) #echarts

#车辆行驶参数模块-误差分析
@app.route('/error',methods=['GET','POST'])
def error():


    return render_template("error_echarts.html") #echarts控件
    # return render_template("error.html",context1=context1,context=context) #plotly画图


#测试评价模块
@app.route('/test_grade',methods=['GET','POST'])

def index():

    return render_template('test_evaluation.html')

@app.route('/grade',methods=['GET','POST'])
def index2():


    return render_template('jieguotest.html')

@app.route('/login', methods=['GET', 'POST'])
def login():


    return render_template('login.html')


# @property
# def password(self):
#     raise "you cant read it"
#     #使用user.password='asda'设置时存入生成的散列密码
# @password.setter
# def password(self, password):
#     self.password_hash = generate_password_hash(password)

# User Register
@app.route('/register', methods=['GET', 'POST'])
def register():


    return render_template('register.html')

# Check if user logged in

def is_logged_in(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        if 'logged_in' in session:
            return f(*args, **kwargs)
        else:
            flash('Unauthorized, Please login', 'danger')
            return redirect(url_for('login'))
    return wrap

# Logout
@app.route('/logout')
@is_logged_in
def logout():

    return redirect(url_for('login'))


if __name__ == '__main__':
#指定服务器IP和端口
    # app=create_app()
    # webbrowser.open('http://127.0.0.1:5000/',new=2,autoraise=False)
    app.run(host='127.0.0.1',port='5000',debug=False,threaded=True)
    # from livereload import Server
    #
    # server = Server(app.wsgi_app)
    # server.watch('**/*.*')
    # server.serve()
