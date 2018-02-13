
switch.condition.menu<-function(){
  switch(letter,

         "S" = {place.size.condition()
           sys.call(which = -1)
           condition.menu()},

         "M" = {place.mig.condition()
           sys.call(which = -1)
           condition.menu()},

         "T" = {place.time.condition()
           sys.call(which = -1)
           condition.menu()},


         "1" = {print(.e$size.matrix)
                print("-----------------")
                sys.call(which = -1)
                condition.menu()},

         "2" = {print(.e$mig.matrix)
                print("-----------------")
                sys.call(which = -1)
                condition.menu()},

         "3" = {print(.e$time.matrix)
                print("-----------------")
                sys.call(which = -1)
                condition.menu()},

         "B" = {sys.call(which = -1)
           main.menu()}

         )}
