import cv2
import os

'''
Change the path in line 15 and execute this programme to
sample 500 of face images in order to make dataset, 
please update folders to the Github after sampling has been done. 
                                                                --Hana
'''

imagecount = 1000
cap = cv2.VideoCapture(0)

#path of classfier should be changed to the corresponding saving path of haarcascade_frontalface_default.xml, this path only exists in Hana's laptop
classfier = cv2.CascadeClassifier("D:\opencv\opencv\sources\data\haarcascades\haarcascade_frontalface_default.xml")

#get current path
current_path = os.getcwd()
path = current_path + '/train5'

def main():
    if not os.path.exists(path):
        os.mkdir(current_path + '/train5')#create a folder called 'train'
    count = 0
    while True:
        ret, frame = cap.read()
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        faceRects = classfier.detectMultiScale(gray, scaleFactor=1.2, minNeighbors=15, minSize=(32, 32))
        
        if len(faceRects) > 0:
            for faceRect in faceRects:
                x, y, w, h = faceRect
                cv2.rectangle(frame, (x - 10, y - 10), (x + w + 10, y + h + 10), (0, 255, 0), 2)
                count += 1
                face_gray = gray[(y - 10):(y + h + 10), (x - 10): (x + w + 10)]
                cv2.imwrite(path + '/' + str(count) + '.jpg' , face_gray)
                print("image {} is stored successfully!".format(count))
        
        cv2.imshow('frame', frame)
        if cv2.waitKey(1) & 0xFF == ord('q') or count >= imagecount:
            break

    cap.release()
    cv2.destroyAllWindows()

if __name__ == '__main__':
    main()