import gc
import abc
from multiprocessing import Process, Manager

import cv2


# 定义抽象基类，此类不能直接实例化
# 做好框架
# 其子类只用实现.process_image方法，返回任意图像算法处理后的从缓存栈中读取的图片
class ABVideoCapture(abc.ABC):
    def __init__(self, cam, top=100):
        self.stack = Manager().list()
        self.max_cache = top
        self.write_process = Process(target=self.__class__.write, args=(self.stack, cam, top))
        self.write_process.start()
        self.__read_gen = self.read_gen()

    @abc.abstractmethod
    def process_image(self, image):
        """
        对输入的图片进行处理并返回处理后的图片
        """
        return self

    def read_gen(self):
        while True:
            if len(self.stack) != 0:
                img = self.process_image(self.stack.pop())
                yield img

    def read(self):
        try:
            return True, next(self.__read_gen)
        except StopIteration:
            return False, None
        except TypeError:
            raise TypeError('{}.read_gen必须为生成器函数'.format(self.__class__.__name__))

    def __iter__(self):
        yield from self.__read_gen

    def release(self):
        self.write_process.terminate()

    def __del__(self):
        self.release()

    @staticmethod
    def write(stack, cam, top):
        """向共享缓冲栈中写入数据"""
        cap = cv2.VideoCapture(cam)
        while True:
            _, img = cap.read()
            if _:
                stack.append(img)
                # 每到一定容量清空一次缓冲栈
                # 利用gc库，手动清理内存垃圾，防止内存溢出
                if len(stack) >= top:
                    del stack[:]
                    gc.collect()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.release()


# 继承ABVideoCapture，对缓存栈中的图片不做处理直接返回
class VideoCapture(ABVideoCapture):
    def process_image(self, image):
        # 这里对图像的处理算法可以随意制定
        # 因为opencv读取的图片并非jpeg格式，因此要用motion JPEG模式需要先将图片转码成jpg格式图片
        ret, jpeg = cv2.imencode('.jpg', image)
        return jpeg.tobytes()
        # return image