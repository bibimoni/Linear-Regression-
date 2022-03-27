import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;

public class InputReader {
	
	File file;
	
	
	public InputReader(File file) {
		this.file = file;
	}
	
	public String readFile() throws FileNotFoundException {
		String str = "";
		int i = 0;
		Scanner sc = new Scanner(file);		
		try {
			while(sc.hasNextLine()) {
				str += " "+sc.nextLine();
				i++;
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}		
		return str;		
	}
	public int getLength() {
		int i = 0;
		try {
			Scanner sc = new Scanner(file);		
			while(sc.hasNextLine()) {
				sc.nextLine();
				i++;
			}
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		
		return i;
	}
}
