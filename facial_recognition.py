import cv2
from convolutional_neural_network import cnn_output


class Person:
    def __init__(self, name, id, health_status) -> None:
        self.name = name
        self.health_status = health_status
        self.id = id

Hana = Person('Hana', '23', 'Bad')
Ben = Person('Ben', '23', 'Good')
Ayoan = Person('Ayoan', '23', 'Good')
Hammond = Person('Hammond', '23', 'Bad')



def facial_recognition():
    cap = cv2.VideoCapture(0)
    classfier = cv2.CascadeClassifier("D:\opencv\opencv\sources\data\haarcascades\haarcascade_frontalface_default.xml")
    while True:
        ret, frame = cap.read()
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        health = 1
        faceRects = classfier.detectMultiScale(frame, scaleFactor=1.2, minNeighbors=10, minSize=(32, 32))
        if len(faceRects) > 0:
            for faceRect in faceRects:
                x, y, w, h = faceRect
                
                face_gray = gray[(y - 10):(y + h + 10), (x - 10): (x + w + 10)]
                #print(size(face_gray))
                probs, index = cnn_output(face_gray)
                print(probs)
                if index == 0:
                    cv2.putText(frame, '{} {} {}'.format(Hana.name, Hana.id, Hana.health_status), (x-15, y-30), cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 0, 255), 2)
                    health = 0
                elif index == 1:
                    cv2.putText(frame, '{} {} {}'.format(Hammond.name, Hammond.id, Hammond.health_status), (x-15, y-30), cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 0, 255), 2)
                    health = 0
                elif index == 2:
                    cv2.putText(frame, '{} {} {}'.format(Ayoan.name, Ayoan.id, Ayoan.health_status), (x-15, y-30), cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 0, 255), 2)
                    health = 1
                elif index == 3:
                    cv2.putText(frame, '{} {} {}'.format(Ben.name, Ben.id, Ben.health_status), (x-15, y-30), cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 0, 255), 2)
                    health = 1
                if health:
                    cv2.rectangle(frame, (x - 10, y - 10), (x + w + 10, y + h + 10), (0, 255, 0), 2)
                else:
                    cv2.rectangle(frame, (x - 10, y - 10), (x + w + 10, y + h + 10), (0, 0, 255), 2)
        cv2.imshow('frame', frame)
        
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break
    cap.release()
    cv2.destroyAllWindows()

if __name__ == '__main__':
    facial_recognition()