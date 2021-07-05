

# import os
#
# from flask import Flask, render_template, Response, make_response
# from camera import VideoCamera
# import numpy as np
# import cv2
# import threading
# from copy import deepcopy
#
# thread_lock = threading.Lock()
# thread_exit = False
# app = Flask(__name__)
# class myThread(threading.Thread):
#     def __init__(self, camera_id, img_height, img_width):
#         super(myThread, self).__init__()
#         self.camera_id = camera_id
#         self.img_height = img_height
#         self.img_width = img_width
#         self.frame = np.zeros((img_height, img_width, 3), dtype=np.uint8)
#
#     def get_frame(self):
#         return deepcopy(self.frame)
#
#     def run(self):
#         global thread_exit
#         cap = cv2.VideoCapture(self.camera_id)
#         while not thread_exit:
#             ret, frame = cap.read()
#             if ret:
#                 frame = cv2.resize(frame, (self.img_width, self.img_height))
#                 thread_lock.acquire()
#                 self.frame = frame
#                 thread_lock.release()
#             else:
#                 thread_exit = True
#         cap.release()
#
# # 相机推流
# def gen():
#     global thread_exit
#     camera_id = "rtsp://admin:kcy123456@172.20.20.102:554/h264/ch1/main/av_stream"
#     img_height = 480
#     img_width = 640
#     thread = myThread(camera_id, img_height, img_width)
#     thread.start()
#
#     while not thread_exit:
#         thread_lock.acquire()
#         frame = thread.get_frame()
#         yield (frame )
#         thread_lock.release()
#     #
#     # thread.join()
#     #
#
#
# # 相机喂流
# @app.route('/video_feed')
# def video_feed():
#     return Response(gen(),
#                     mimetype='multipart/x-mixed-replace; boundary=frame')
#
#
# # 当前实时相机画面
# @app.route('/')
# def cur_camera():
#     return render_template('camera.html')
#
#
# if __name__ == '__main__':
#     app.run(host='127.0.0.1',port="8000", debug=True)


import os
import cv2
from flask import Flask, render_template, Response, make_response
from camera import VideoCamera



app = Flask(__name__)


# 相机推流
def gen(camera):
    while True:
        frame = camera.read()
        yield (b'--frame\r\n'
               b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n\r\n')


# 相机喂流
@app.route('/video_feed')
def video_feed():
    return Response(gen(VideoCapture("rtsp://admin:kcy123456@172.20.20.102:554/h264/ch1/main/av_stream")),
                    mimetype='multipart/x-mixed-replace; boundary=frame')


# 当前实时相机画面
@app.route('/')
def cur_camera():
    return render_template('camera.html')


if __name__ == '__main__':
    app.run(host='127.0.0.1',port="8000", debug=True)